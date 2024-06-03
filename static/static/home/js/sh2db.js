function showAlert(message, alerttype) {
  // Alerttype: success, info, warning, danger
  // See https://getbootstrap.com/docs/3.3/components/#alerts
  $("#sh2db_alert_placeholder").append('<div id="sh2db_alert" class="alert alert-' +
    alerttype + '"><button type="button" class="close" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
    message + "</div>");

  setTimeout(function() {
    $("#sh2db_alert").remove();
  }, 4000);
}