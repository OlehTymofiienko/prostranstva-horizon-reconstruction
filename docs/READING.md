# Языки и темы чтения / Languages and reading themes

Дополнение от 23 сентября 2026 года к научному выпуску 1.0. Последний научный этап — HF4-B5. Подтверждённые расчёты не повторялись.

В верхней панели каждой страницы расположены ссылки «Русский / English» и выбор темы: системная, светлая, тёмная или тёплая. Тема сохраняется в браузере. Ссылки на английских страницах ведут к соответствующим английским материалам; исходные численные файлы и программы доступны по прежним адресам.

Английская версия включает главную страницу, каталог, четыре интерактивных просмотра, 18 отчётов и README. Отчёты обозначены как English reading edition: это редакционные английские версии с основными методами, результатами и ограничениями; исходный русский отчёт доступен рядом. Подписи внутри ранее созданных изображений и видео, архивные программы и названия исходных файлов сохранены. Научные цветовые шкалы и данные не изменяются при выборе темы.

## Rebuilding the reading layer

The English reading editions are in `docs/en/`; interface translations are in `scripts/reading-ui.json` and `scripts/reading-dynamic.json`. The renderer uses the vendored Markdown parser. No scientific evolution or fitting is run by these commands:

```bash
python scripts/localize_site.py
python scripts/check_release.py
python scripts/check_reading.py
```

`dist/` is the deployed static site. The existing manual GitHub Pages workflow publishes this directory with `contents: read`. No additional workflow permissions, automated commits or releases are introduced.

The English editions preserve the reported results and limitations, with historical pending actions explicitly identified as historical. The final project report takes precedence over early-stage parameter interpretations. All 18 Russian source reports remain available. Existing figure/video text, numerical archives and program comments are retained in their original language.

Theme selection affects reading surfaces and interface contrast. It does not recolour scientific field maps, alter arrays, change thresholds or rerun numerical analysis. Dark and warm modes are reading preferences, not medical claims.

Validation records are in `docs/reading-validation.json` and `docs/release-validation.json`. Publication status and continuation state are recorded in `docs/PUBLICATION_STATUS.md` and `CHECKPOINT.json`.
