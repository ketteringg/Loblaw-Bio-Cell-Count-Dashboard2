CREATE TABLE subjects (
    subject_id   TEXT PRIMARY KEY,
    project      TEXT NOT NULL,
    condition    TEXT NOT NULL,
    age          INTEGER NOT NULL,
    sex          TEXT NOT NULL,
    treatment    TEXT NOT NULL,
    response     TEXT
);

CREATE TABLE samples (
    sample_id                  TEXT PRIMARY KEY,
    subject_id                 TEXT NOT NULL REFERENCES subjects(subject_id),
    sample_type                TEXT NOT NULL,
    time_from_treatment_start  INTEGER NOT NULL
);

CREATE TABLE cell_counts (
    sample_id   TEXT NOT NULL REFERENCES samples(sample_id),
    population  TEXT NOT NULL,
    count       INTEGER NOT NULL,
    PRIMARY KEY (sample_id, population)
);
