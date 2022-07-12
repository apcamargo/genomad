# Copyright (c) 2020 ReDNA Labs Co., Ltd.

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import tensorflow as tf

from tensorflow.keras.layers import (
    Conv1D,
    LeakyReLU,
    SpatialDropout1D,
    Concatenate,
    Dropout,
    Layer,
)

from tensorflow.keras.regularizers import l2

import sys
import numpy as np


def IGLOO1D_Block(
    incoming_layer,
    nb_patches,
    nb_filters_conv1d,
    patch_size=4,
    padding_style="causal",
    nb_stacks=1,
    conv1d_kernel=3,
    dropout_rate=0.01,
    l2_reg=0.000001,
    transformer_style=False,
    spatial_dropout=True,
    pooling_size=1,
    incoming_proj=0,
):

    x = Conv1D(nb_filters_conv1d, conv1d_kernel, padding=padding_style)(incoming_layer)
    x = LeakyReLU(alpha=0.1)(x)
    x = (
        SpatialDropout1D(dropout_rate)(x)
        if spatial_dropout
        else Dropout(dropout_rate)(x)
    )
    x_igloo = IGLOO1D_kernel(
        patch_size,
        nb_patches,
        dropout_rate,
        l2_reg,
        transformer_style,
        pooling_size,
        incoming_proj,
    )(x)
    layers = [x_igloo]
    if nb_stacks > 1:
        for _ in range(nb_stacks - 1):
            x = Conv1D(nb_filters_conv1d, conv1d_kernel, padding=padding_style)(x)
            x = LeakyReLU(alpha=0.1)(x)
            x = (
                SpatialDropout1D(dropout_rate)(x)
                if spatial_dropout
                else Dropout(dropout_rate)(x)
            )
        x_igloo = IGLOO1D_kernel(
            patch_size,
            nb_patches,
            dropout_rate,
            l2_reg,
            transformer_style,
            pooling_size,
            incoming_proj,
        )(x)
        layers.append(x_igloo)
    return Concatenate()(layers) if nb_stacks > 1 else layers[0]


class IGLOO1D_kernel(Layer):
    def __init__(
        self,
        patch_size,
        nb_patches,
        dropout_rate,
        l2_reg,
        transformer_style,
        pooling_size,
        incoming_proj,
    ):
        super(IGLOO1D_kernel, self).__init__()
        self.patch_size = patch_size
        self.nb_patches = nb_patches
        self.dropout_rate = dropout_rate
        self.l2_reg = l2_reg
        self.pooling_size = pooling_size
        self.incoming_proj_dim = incoming_proj
        self.transformer_style = transformer_style

    def core_igloo_layer_initializer(self, shape, dtype=None):
        M = gen_filters_igloo(
            self.patch_size,
            self.nb_patches,
            self.vector_size,
            build_backbone=False,
            return_sequences=False,
        )
        M.astype(int)
        return M

    def build(self, input_shape):
        self.batch_size = input_shape[0]
        self.vector_size = input_shape[1]
        if self.incoming_proj_dim > 0:
            self.w_incoming = self.add_weight(
                shape=(input_shape[2], self.incoming_proj_dim),
                initializer="glorot_uniform",
                trainable=True,
                regularizer=l2(self.l2_reg),
                name="w_incoming",
            )
        self.num_channels_input = input_shape[2]
        self.patches = self.add_weight(
            shape=(int(self.nb_patches), self.patch_size, 1),
            initializer=self.core_igloo_layer_initializer,
            trainable=False,
            name="random_patches",
            dtype=np.int32,
        )
        if self.incoming_proj_dim == 0:
            self.w_mult = self.add_weight(
                shape=(1, self.nb_patches, self.patch_size, self.num_channels_input),
                initializer="glorot_uniform",
                trainable=True,
                regularizer=l2(self.l2_reg),
                name="w_mult",
            )
            self.w_summer = self.add_weight(
                shape=(1, self.patch_size * self.num_channels_input, 1),
                initializer="glorot_uniform",
                trainable=True,
                regularizer=l2(self.l2_reg),
                name="w_summer",
            )
        else:
            self.w_mult = self.add_weight(
                shape=(1, self.nb_patches, self.patch_size, self.incoming_proj_dim),
                initializer="glorot_uniform",
                trainable=True,
                regularizer=l2(self.l2_reg),
                name="w_mult",
            )
            self.w_summer = self.add_weight(
                shape=(1, self.patch_size * self.incoming_proj_dim, 1),
                initializer="glorot_uniform",
                trainable=True,
                regularizer=l2(self.l2_reg),
                name="w_summer",
            )
        self.w_bias = self.add_weight(
            shape=(1, self.nb_patches),
            initializer="glorot_uniform",
            trainable=True,
            regularizer=l2(self.l2_reg),
            name="w_bias",
        )
        if self.transformer_style:
            self.w_qk = self.add_weight(
                shape=(self.nb_patches, int(self.vector_size / self.pooling_size)),
                initializer="glorot_uniform",
                trainable=True,
                regularizer=l2(self.l2_reg),
                name="w_qk",
            )

            self.w_v = self.add_weight(
                shape=(1, self.num_channels_input, self.num_channels_input),
                initializer="glorot_uniform",
                trainable=True,
                regularizer=l2(self.l2_reg),
                name="w_v",
            )

    def call(self, y):
        y_next = tf.matmul(y, self.w_incoming) if self.incoming_proj_dim > 0 else y
        M = tf.transpose(y_next, [1, 2, 0])
        M = tf.gather_nd(M, self.patches)
        mpi = tf.transpose(M, [3, 0, 1, 2])
        mpi = tf.multiply(self.w_mult, mpi)
        if self.incoming_proj_dim == 0:
            mpi = tf.reshape(
                mpi, [-1, self.nb_patches, self.patch_size * self.num_channels_input]
            )
        else:
            mpi = tf.reshape(
                mpi, [-1, self.nb_patches, self.patch_size * self.incoming_proj_dim]
            )
        mpi = tf.matmul(mpi, self.w_summer)
        mpi = tf.squeeze(mpi, axis=-1)
        mpi = mpi + self.w_bias
        if self.transformer_style:
            y_proj = tf.matmul(y, self.w_v)
            if self.pooling_size > 1:
                y_proj = tf.keras.layers.MaxPool1D(pool_size=self.pooling_size)(y_proj)
            alpha = tf.matmul(mpi, self.w_qk)
            alpha = tf.nn.softmax(alpha)
            mpi = tf.matmul(tf.expand_dims(alpha, axis=1), y_proj)
            mpi = tf.squeeze(mpi, axis=1)
        else:
            mpi = LeakyReLU(alpha=0.1)(mpi)
        return mpi


def gen_filters_igloo(
    patch_size,
    nb_patches,
    vector_size,
    return_sequences,
    build_backbone=True,
    consecutive=False,
    nb_sequences=-1,
):
    outa = []
    vector_size = int(vector_size)
    for step in range(vector_size):
        if (step != vector_size - 1) and (return_sequences == False):
            continue
        if (
            return_sequences == True
            and (nb_sequences != -1)
            and step < vector_size - nb_sequences
        ):
            continue
        collect = []
        if step < patch_size:
            for _ in range(nb_patches):
                randy_H = np.random.choice(range(step + 1), patch_size, replace=True)
                first = [[pp] for pp in randy_H]
                collect.append(first)
        else:
            sorting = True
            if build_backbone:
                maximum_its = int((step / (patch_size - 1)) + 1)
                if maximum_its > nb_patches:
                    print("nb_patches too small, recommended above:", maximum_its)
                    sys.exit()
                for jj in range(maximum_its):
                    if iter == 0:
                        randy_H = [step - pp for pp in range(patch_size)]
                    else:
                        randy_H = [
                            max(step - (jj * (patch_size - 1)) - pp, 0)
                            for pp in range(patch_size)
                        ]
                    first = [[pp] for pp in randy_H]
                    collect.append(first)
                rest_iters = max(nb_patches - maximum_its, 0)
                for _ in range(rest_iters):
                    if not consecutive:
                        randy_B = np.random.choice(
                            range(step + 1), patch_size, replace=False
                        )
                    else:
                        uniq = np.random.choice(
                            range(max(0, step + 1 - patch_size + 1)),
                            1,
                            replace=False,
                        )
                        randy_B = [uniq[0] + pp for pp in range(patch_size)]
                    if sorting:
                        randy_B = sorted(randy_B)
                    first = [[pp] for pp in randy_B]
                    collect.append(first)
            else:
                for _ in range(nb_patches):
                    if not consecutive:
                        randy_B = np.random.choice(
                            range(step + 1), patch_size, replace=False
                        )
                    else:
                        uniq = np.random.choice(
                            range(max(0, step + 1 - patch_size + 1)),
                            1,
                            replace=False,
                        )
                        randy_B = [uniq[0] + pp for pp in range(patch_size)]
                    if sorting:
                        randy_B = sorted(randy_B)
                    first = [[pp] for pp in randy_B]
                    collect.append(first)
            collect = np.stack(collect)
        outa.append(collect)
    outa = np.stack(outa)
    if return_sequences == False:
        outa = np.squeeze(outa, axis=0)
    return outa


class BranchAttention(Layer):
    def __init__(self):
        super(BranchAttention, self).__init__()

    def build(self, input_shape):
        self.batch_size = input_shape[0][0]
        self.fulloutput = input_shape[2][1]
        self.w_1 = self.add_weight(
            shape=(1, 6),
            initializer="glorot_normal",
            trainable=True,
            name="w_1",
        )
        self.w_2 = self.add_weight(
            shape=(1, 6), initializer="glorot_normal", trainable=True, name="w_2"
        )
        super(BranchAttention, self).build(input_shape)

    def call(self, inputs):
        weighting_parameter = inputs[0]
        first_branch = inputs[1]
        second_branch = inputs[2]
        alpha = tf.matmul(weighting_parameter, self.w_1) + self.w_2
        weighted_first_branch = tf.multiply(alpha[:, 0:3], first_branch)
        weighted_second_branch = tf.multiply(alpha[:, 3:6], second_branch)
        return tf.reduce_mean([weighted_first_branch, weighted_second_branch], axis=0)

    def compute_output_shape(self, input_shape):
        return input_shape[0][0], self.fulloutput
