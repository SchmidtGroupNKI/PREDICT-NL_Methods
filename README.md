# Improving the predictive performance of Cox models: A framework for updating and extension with new markers

Continuous monitoring and improvement of clinical prediction models are essential for maintaining their performance. We developed a systematic framework for updating and extension of a Cox regression model for time-to-event data, with illustration for the PREDICT model. PREDICT is a widely-used model for deciding the need for adjuvant therapy for breast cancer. Like any prediction model, it may underperform in new patients. We illustrated the updating and extension of PREDICT for the MINDACT trial cohort, which comprised somewhat more selectively recruited participants treated according to more recent clinical guidelines compared to the registry data underlying the development of PREDICT.

# Syntax files
| File                   | Description             |
| :----                  | :----                   |
| 1_Cleaning.Rmd                       | Cleaning, preparation, and imputation of the MINDACT dataset.
| 2_Descriptive.Rmd                    | Descriptive tables.
| 3_PREDICT_validation.Rmd             | Validation of the PREDICT v2.3 model.
| 4_PREDICT_updating.Rmd               | PREDICT extension and validation of the extended models.
| 4_PREDICT_updating_continued.Rmd     | PREDICT extension and validation of the extended models (calibration plots, decision curves, and net benefit).
| 5_Method_selection.Rmd               | Application of the closed test procedure on the full MINDACT dataset.
| 5_Method_selection_simulation.Rmd    | Simulation of the closed test procedure on samples of the MINDACT dataset.

# Contact
Mary Ann E. Binuya <br/>
Netherlands Cancer Institute <br/>
[m.binuya@nki.nl](m.binuya@nki.nl)

# Authors
| Author                 | Role   | Description             |
| :----                  | :----: | :----                   |
| Mary Ann Binuya   | Author | Development and support |
| Ellen Engelhardt  | Author | Review                  |
| Terry Chan        | Author | Development and review  |
| Martijn Heymans   | Author | Review                  |
| Marjanka Schmidt  | Author | Review                  |
| Ewout Steyerberg  | Author | Development and review  |
