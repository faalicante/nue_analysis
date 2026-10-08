# CNN 2+1D a due teste

`training/train_cnn21d.py` implementa la CNN fattorizzata, la loss multi-task, checkpoint,
early stopping, metriche e diagnostiche. Il task iniziale deriva esclusivamente
da `sample_type`: `2` identifica il signal, `1` l'hard e `0` il Poisson.
Nella configurazione corrente il signal è positivo per la BCE solo sopra 10 mrad,
ma tutti i signal (incluso `signal_l5.root`) sono usati dalla regressione. Hard e
una quota del 10% di Poisson sono negativi per la BCE. L'input è `[B,1,z,y,x]`; `slope_xy`
segue gli assi ROOT verificati in `CNN_DATALOADER.md`, quindi la direzione è
`phi = atan2(sy,sx)`.

Esecuzioni sicure e separate:

```sh
/tmp/crop-bw-venv/bin/python -m pytest -q visualization/tests/test_cnn21d.py
/tmp/crop-bw-venv/bin/python -m training.train_cnn21d --config training/configs/cnn21d.yaml --mode tiny
/tmp/crop-bw-venv/bin/python -m training.train_cnn21d --config training/configs/cnn21d.yaml --mode pilot
```

Il training completo non viene avviato implicitamente. Quando autorizzato:

```sh
/tmp/crop-bw-venv/bin/python -m training.train_cnn21d --config training/configs/cnn21d.yaml --mode full
```

Gli output sono separati sotto `runs/cnn21d/{tiny,pilot,full}`. La risoluzione
angolare riportata è la deviazione standard (popolazione) del residuo
`theta_pred-theta_true`; l'errore azimutale usa il residuo circolare in `[-pi,pi]`.
La diagnostica edge/central usa, in assenza di una coordinata di impatto nei
metadati, i quartili della frazione di attivazione positiva contenuta nel bordo
xy largo 3 pixel.

## Training bilanciato e checkpoint separati

`training/configs/cnn21d_signal_p24_l5_dual_checkpoint_pilot.yaml` abilita:

- regressione pesata per bin di theta con pesi inverse-square-root limitati;
- layer-drop train-only (12% un layer, 3% due layer consecutivi);
- quota train `Poisson:hard:signal = 10:40:50` nei subset;
- `best_classification_model.pt`, selezionato sulla AUPRC;
- `best_regression_model.pt`, selezionato sul MAE angolare medio tra bin.

Nello scan i due checkpoint vanno passati insieme:

```sh
python3 -m scanning.scan_cnn21d_volumes \
  --input INPUT.h5 \
  --config training/configs/cnn21d_signal_p24_l5_dual_checkpoint_pilot.yaml \
  --checkpoint RUN/best_classification_model.pt \
  --regression-checkpoint RUN/best_regression_model.pt \
  --output OUTPUT.h5 --crop-size 20 --stride 10
```

Lo score viene dal checkpoint di classificazione; `slope_x`, `slope_y`, theta e
phi vengono dal checkpoint di regressione. `best_model.pt` resta un alias del
checkpoint di classificazione per compatibilità con i comandi precedenti.
