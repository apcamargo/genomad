import tensorflow as tf

from keras import layers as kl
from keras import Model, Layer

from genomad.neural_network import igloo


class OneHotLayer(Layer):
    def call(self, x):
        return tf.one_hot(x, depth=257, axis=-1)


def create_encoder():
    inputs = kl.Input(shape=(5_997,), dtype="int64")
    embedded = OneHotLayer()(inputs)
    outputs = igloo.IGLOO1D_Block(
        embedded,
        nb_patches=2100,
        nb_filters_conv1d=128,
        nb_stacks=3,
        conv1d_kernel=6,
        pooling_size=8,
        dropout_rate=0.2,
        l2_reg=1e-3,
        transformer_style=True,
    )
    outputs = kl.Dense(512)(outputs)
    outputs = kl.BatchNormalization()(outputs)
    outputs = kl.Activation("relu")(outputs)
    return Model(inputs=inputs, outputs=outputs)


def create_classifier():
    encoder = create_encoder()
    for layer in encoder.layers:
        layer.trainable = False
    inputs = kl.Input(shape=(5_997,))
    features = encoder(inputs)
    features = kl.Dense(512)(features)
    features = kl.BatchNormalization()(features)
    features = kl.Activation("relu")(features)
    features = kl.Dropout(0.2)(features)
    outputs = kl.Dense(3, activation="softmax")(features)
    return Model(inputs=inputs, outputs=outputs)
