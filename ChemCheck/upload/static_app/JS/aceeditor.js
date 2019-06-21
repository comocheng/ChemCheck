$(document).ready(function(){
    $("#edit").click(function(){
      $.get("{{ mech.ck_mechanism_file.url }}");
    });

    function initAceEditor() {
        var editor = ace.edit("aceEditor");
        editor.setTheme("ace/theme/solarized_light");
        editor.getSession().setMode("ace/mode/html");
        editor.setFontSize("15px") ;
        editor.setPrintMarginColumn(false);
        editor.setValue($("#file-content").text());
        editor.setOptions({
        maxLines: 2000
        });