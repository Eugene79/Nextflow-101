function onRender(input) {

    var validationErrors = [];
    var validationWarnings = [];

    if (input.currentAnalysisSettings === null) {
        //null first time, to use it in the remainder of he javascript
        input.currentAnalysisSettings = input.analysisSettings;
    }

    switch(input.context) {
        case 'Initial': {
            renderInitial(input, validationErrors, validationWarnings);
            break;
        }
        case 'FieldChanged': {
            renderFieldChanged(input, validationErrors, validationWarnings);
            break;
        }
        case 'Edited': {
            renderEdited(input, validationErrors, validationWarnings);
            break;
        }
        default:
            return {};
    }

    return {
        'analysisSettings': input.currentAnalysisSettings,
        'settingValues': input.settingValues,
        'validationErrors': validationErrors,
        'validationWarnings': validationWarnings
    };
}

function renderInitial(input, validationErrors, validationWarnings) {
}

function renderEdited(input, validationErrors, validationWarnings) {
}

function renderFieldChanged(input, validationErrors, validationWarnings) {
}

function findField(input, fieldId){
    var fields = input.currentAnalysisSettings['fields'];
    for (var i = 0; i < fields.length; i++){
        if (fields[i].id === fieldId) {
            return fields[i];
        }
    }
    return null;
}
