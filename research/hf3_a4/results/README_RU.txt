BH-YITIAN-A4

Начните с BH_YITIAN_A4_report.md.
Численные JSON содержат все варианты, а не только успешные.
BH_HF3_resume_state_A4.json — точка продолжения; A1–A3 и эволюция не повторяются.

Повторение с исходным архивом Yitian:
python BH_YITIAN_A4.py --input /путь/common_horizon_moments_data_and_scripts.zip --output repeated

Частоты уже включены. Их пересчёт: добавить --recompute-spectrum (нужен qnm==0.4.4).
Только оформление готовых JSON: --render-only с прежней папкой --output.
Эволюция Einstein Toolkit/SpEC не запускается.
