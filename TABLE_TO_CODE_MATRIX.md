# Table-to-code matrix

Exact commands and paths: README.md. Fixed configurations: specs/ and inputs/protocol/. Full replay outputs were compared with frozen references, then manuscript rows.

|Target|Status|Route|
|---|---|---|
|Main1|PASS|experiment.py --dataset mechanism → summarize.py → build_tables.py|
|Main2|PASS|experiment.py --dataset lucas → summarize.py → build_tables.py|
|Main3|PASS|experiment.py --dataset lucas → summarize.py → build_tables.py|
|S1|PASS|experiment.py --dataset lucas → summarize.py → build_tables.py|
|S2|PASS|experiment.py --dataset lucas → summarize.py → build_tables.py|
|S3|PASS|experiment.py --dataset lucas → summarize.py → build_tables.py|
|S4|PASS|gap_run.py → gap_summarize.py → check_tables.py|
|S5|PASS|directional.py --fields 128 --labels 128 --with-null → check_tables.py|
|S6|PASS|directional.py --fields 128 --labels 128 --with-null → check_tables.py|
|S7|PASS|experiment.py --dataset lucas → crosscondition.py → check_tables.py|
|SVD/derivative/delete/remainder|PASS|audit_regression.py|
|theory/code identities|PASS|theory_bridge.py|
|slope boundary|PASS|slope_boundary.py|
|source frame|PASS|rebuild_frame.py → compare_additional.py|
