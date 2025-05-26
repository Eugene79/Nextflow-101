function onSubmit(input) {
    var validationErrors = [];

    return {
        'settings': input.settings,
        'validationErrors': validationErrors
    };
}
