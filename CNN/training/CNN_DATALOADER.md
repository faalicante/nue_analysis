# PyTorch Dataset e DataLoader per CNN 2+1D

## Convenzione verificata degli assi

La convenzione non è dedotta da `imshow`. Nel produttore `crop_bw.C`, per ogni
slice z vengono letti i bin ROOT con `GetBinContent(sourceX, sourceY)` e scritti
nell'array C all'indice lineare

```text
(z * 32 * 32) + (y * 32) + x
```

Il convertitore non traspone l'array e salva
`regions_raw[event,z,y,x]`. Quindi in NumPy:

- `volume[z,y,x]`: z è la slice, y è la riga, x è la colonna;
- `slope_xy=(slope_x,slope_y)` segue gli assi fisici ROOT x e y;
- le figure usano `origin="lower"`, perciò x cresce a destra e y cresce in alto.

Con questa combinazione, una rotazione fisica antioraria di 90 gradi non è
`np.rot90(k=+1)`: è `np.rot90(k=-1, axes=(-2,-1))`. Un punto sintetico posto a
destra del centro (+x) finisce sopra il centro (+y), e la slope segue
`(sx,sy) -> (-sy,sx)`. I test verificano sia il punto sia il vettore.

Le trasformazioni sono applicate in quest'ordine:

1. rotazione fisica xy antioraria di `k * 90` gradi;
2. `reflect_x`: `x -> -x`, inversione delle colonne, `sx -> -sx`;
3. `reflect_y`: `y -> -y`, inversione delle righe, `sy -> -sy`.

## Pipeline

`HDF5CNN21DDataset` usa solo `regions_raw`, cioè conteggi interi grezzi, non
smussati. Per il training estrae mediante slicing un crop comune a tutte le
slice z. Con il default `max_jitter=5`, il crop parte da `(6+dx,6+dy)` con
`dx,dy` interi uniformi in `[-5,+5]`. Non sono usati `roll`, padding o wrapping.

`background_mu` nasce nel TTree ROOT come scalare per evento. Il convertitore
corrente lo serializza ripetendolo sulle 57 posizioni z; il Dataset verifica che
le 57 copie siano identiche e lo riconduce a uno scalare. Subito dopo il crop
viene applicata

```text
X[z,y,x] = (H[z,y,x] - background_mu) / sqrt(background_mu + 1e-6)
```

Non viene applicata alcuna soglia a 5 sigma. Rotazioni e riflessioni seguono la
normalizzazione. Validation e test hanno sempre `dx=dy=k=0` e nessuna
riflessione.

Ogni elemento è un dizionario con:

- `volume`: `torch.float32 [1,57,20,20]`;
- `sample_type`: `torch.int64 []`, target a tre classi: 0 Poisson, 1 hard,
  2 signal;
- `slope_xy`: `torch.float32 [2]`;
- `metadata`: indici HDF5/dataset, split, tipo di campione, identificativi
  evento/tile/crop, provenienza ROOT, parametri di augmentation, processo/worker
  e diagnostica di visibilità. Il vecchio `presence` binario resta qui solo per
  diagnosi.

## Jitter ridotto e audit dei positivi

`positive_max_jitter` può ridurre il jitter solo per `sample_type` 1 o 2, per esempio:

```python
config = CNNDataConfig(max_jitter=5, positive_max_jitter=3)
```

Il target `sample_type` non viene mai cambiato. La diagnostica confronta la massa
`max(raw-background_mu, 0)` mantenuta dal crop traslato con quella del crop
centrale; non usa una soglia sigma. Il flag per elemento è
`positive_may_be_insufficiently_visible`. Per contare i positivi che potrebbero
scendere sotto la frazione configurata in almeno una traslazione consentita:

```python
dataset = HDF5CNN21DDataset("cnn_dataset.h5", split="train", config=config)
print(dataset.audit_positive_visibility())
```

La metrica è volutamente una segnalazione QA, non una ridefinizione della label.

## Uso

Installazione nell'ambiente già usato per il convertitore:

```sh
python -m pip install -r requirements-cnn.txt
```

Creazione dei tre loader:

```python
from training.cnn_dataset import CNNDataConfig, CNNLoaderConfig, create_cnn_dataloaders

loaders = create_cnn_dataloaders(
    "cnn_dataset.h5",
    data_config=CNNDataConfig(max_jitter=5, seed=12345),
    loader_config=CNNLoaderConfig(batch_size=32, num_workers=4, seed=12345),
)
batch = next(iter(loaders["train"]))
```

Il file HDF5 non resta aperto nel costruttore. Ogni processo DataLoader apre
lazy il proprio handle alla prima lettura; handle ereditati da un fork sono
riconosciuti tramite PID/worker e riaperti. Per una nuova augmentation
deterministica a ogni epoca, chiamare `dataset.set_epoch(epoch)` prima di creare
l'iteratore. Il default `persistent_workers=False` rende l'aggiornamento visibile
ai worker a ogni epoca.

Il seed del dataset determina augmentation per `(seed, epoch, hdf5_index)`; il
seed del loader determina shuffle e seed dei worker. A parità di entrambi,
ordine e tensori sono riproducibili.

Le trasformazioni casuali sono abilitate esclusivamente nel Dataset `train`.
Validation e test restituiscono sempre il crop centrale originale. Nel train la
combinazione identità è inclusa nell'estrazione casuale; con i default l'esempio
esattamente originale ha probabilità `1 / (11² * 4 * 2 * 2) = 1/1936` per
lettura. Non viene quindi aggiunta automaticamente una seconda copia originale
di ogni elemento.

## Test e QA visuale

```sh
python -m pytest -q visualization/tests/test_cnn_dataset.py
python -m visualization.smoke_test_cnn_dataloader cnn_dataset.h5 --batch-size 8 --num-workers 2
python -m visualization.qa_cnn_dataset cnn_dataset.h5 --output-dir cnn_dataloader_qa --count 6
```

Lo script seleziona esempi stratificati, mostra la proiezione xy normalizzata
del crop centrale e quella dopo augmentation e sovrappone la slope. La freccia
è lo spostamento corrispondente a 20 bin z (limitato graficamente se esce dal
pannello); i componenti numerici non limitati sono riportati nella figura.
