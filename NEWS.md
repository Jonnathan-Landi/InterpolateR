# Version 1.3-1 (Development)

### New Features

### Bug Fixed

# Version 1.2-0 (Development)

### New Features

-   in the validation mode, the 'mod_validation' parameter has been removed, as it was unnecessary. Now, if the 'stat_validation' parameter is entered, it will validate with the stations entered. If 'stat_validation' is NULL a random validation will be used.
-   We have added the option to manually insert stations for validation, option 'stat_validation', when this parameter is entered it will be used to validate the stations entered by the user, otherwise when stat_validation = NULL, the validation will be performed randomly.
-   We have optimized the RFmerge algorithm and its validation algorithm has been connected to 'validation_module'.
-   We have optimized the RFplus algorithm and its validation algorithm has been connected to 'validation_module'.
-   We have created the 'validation_module' which will contain all the logic for point to pixel validation. All integration methods that require validation will be connected to this module.
-   Add a new interpolation method. RFplus proposed by Landi et al. (2025).
-   Add a new interpolation method. RFmerge proposed by Baez-Villanueva et al. (2020).
-   Add interpolation method Inverse Distance Weighted Interpolation (IDW).

### Bug Fixed

-   An error has been corrected when manually selecting stations for validation, since if the training was not different from 1 the validation was not performed.
-   We have corrected a bug in the validation methods, when the 'random' method was set, the IDs did not match each station in situ. We have made a new approach in the allocation, correcting this problem.
