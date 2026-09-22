# HF4-B5

Этап завершён. Читать `results/BH_HF4_B5_report.md` и
`results/BH_HF4_resume_state_B5.json`. Подтверждённые этапы не пересчитывать.

Для намеренного независимого воспроизведения из каталога `hf4_b5`:

```bash
python -m pip install -r requirements.txt
python BH_HF4_B5.py --input input/BH_HF4_B4_fields_geometry_masks.npz --protocol input/BH_HF4_B5_protocol.json --geometry input/BH_HF4_B4_geometry.csv
python finish_B5.py
```

Нужен Python 3.10+; версии исходного запуска перечислены в definition.json.
Новая эволюция Einstein Toolkit не нужна. Веса площади и исходные маски
читаются из B4. Вектор невязки раскладывается в координатном базисе;
квадрат его нормы не является физической энергией.

Результат B5 описательный. Разрешение симуляции, выбор координат и тетрады
не варьировались. В первом кадре Re Ψ₂ два узла дают основную долю квадрата
невязки; они сохранены и не исключались из расчёта.
