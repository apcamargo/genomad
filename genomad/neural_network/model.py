import tensorflow as tf
from tensorflow.keras import Model
from tensorflow.keras.layers import (
    Activation,
    BatchNormalization,
    Dense,
    Dropout,
    Input,
)
from genomad.neural_network import igloo


def create_encoder():
    inputs = Input(shape=5_997, dtype="int64")
    embedded = tf.one_hot(inputs, depth=257, axis=-1)
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
    outputs = Dense(512)(outputs)
    outputs = BatchNormalization()(outputs)
    outputs = Activation("relu")(outputs)
    return Model(inputs=inputs, outputs=outputs)


def create_classifier():
    encoder = create_encoder()
    for layer in encoder.layers:
        layer.trainable = False
    inputs = Input(shape=5_997)
    features = encoder(inputs)
    features = Dense(512)(features)
    features = BatchNormalization()(features)
    features = Activation("relu")(features)
    features = Dropout(0.2)(features)
    outputs = Dense(3, activation="softmax")(features)
    return Model(inputs=[inputs], outputs=outputs)
