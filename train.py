import json
from pathlib import Path
import numpy as np
import tensorflow as tf

from src.model_loader import ModelLoader
from src.data_generator import DataGenerator  
from preprocessing import Pipeline
from parameters import PARAMETERS

# parameters 
SEED   = PARAMETERS["SEED"]
KMER = PARAMETERS["KMER"]
MODEL  = PARAMETERS["MODEL"]
CLADES = PARAMETERS["CLADES"]
BATCH_SIZE = PARAMETERS["BATCH_SIZE"]
EPOCHS = PARAMETERS["EPOCHS"]
WEIGHTS_PATH = PARAMETERS["WEIGHTS_PATH"]
PREPROCESSING = PARAMETERS["PREPROCESSING"]

# set seed for reproducibility
tf.random.set_seed(SEED)
np.random.seed(SEED)

# -1- Model selection
loader = ModelLoader()
model  = loader(
            model_name=MODEL,
            n_outputs=len(CLADES),
            weights_path=WEIGHTS_PATH
            ) # get compiled model from ./supervised_dna/models

model.compile(optimizer=tf.keras.optimizers.Adam(),
            loss="categorical_crossentropy",
            metrics=["accuracy"]
)

preprocessing = Pipeline(PREPROCESSING)
Path("data/train").mkdir(exist_ok=True, parents=True)
preprocessing.asJSON("data/train/preprocessing.json")

# -2- Datasets
# load list of images for train and validation sets
with open("datasets.json","r") as f:
    datasets = json.load(f)
list_train = datasets["train"]
list_val   = datasets["val"]

## prepare datasets to feed the model
# Instantiate DataGenerator for training set
ds_train = DataGenerator(
    list_train,
    order_output_model = CLADES,
    batch_size = BATCH_SIZE,
    shuffle = True,
    kmer = KMER,
    preprocessing = preprocessing,
)

# Instantiate DataGenerator for validation set
ds_val = DataGenerator(
    list_val,
    order_output_model = CLADES,
    batch_size = BATCH_SIZE,
    shuffle = False,
    kmer = KMER,
    preprocessing = preprocessing,
) 

# -3- Training
# - Callbacks
# checkpoint: save best weights
Path("data/train/checkpoints").mkdir(exist_ok=True, parents=True)
cb_checkpoint = tf.keras.callbacks.ModelCheckpoint(
    filepath='data/train/checkpoints/model-{epoch:02d}-{val_accuracy:.3f}.hdf5',
    monitor='val_loss',
    mode='min',
    save_best_only=True,
    verbose=1
)

# reduce learning rate
cb_reducelr = tf.keras.callbacks.ReduceLROnPlateau(
    monitor='val_loss',
    mode='min',
    factor=0.1,
    patience=8,
    verbose=1,
    min_lr=0.00001
)

# stop training if
cb_earlystop = tf.keras.callbacks.EarlyStopping(
    monitor='val_loss',
    mode='min',
    min_delta=0.001,
    patience=10,
    verbose=1
)

# save history of training
Path("data/train").mkdir(exist_ok=True, parents=True)
cb_csvlogger = tf.keras.callbacks.CSVLogger(
    filename='data/train/training_log.csv',
    separator=',',
    append=False
)

model.fit(
    ds_train,
    validation_data=ds_val,
    epochs=EPOCHS,
    callbacks=[
        cb_checkpoint,
        cb_reducelr,
        cb_earlystop,
        cb_csvlogger,
        ]
)