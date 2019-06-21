$(document).ready(function(){
    $(".btn btn-primary btn-sm").click(function(e){
      var contents = e.target.result;
      displayContents(contents)
      function displayContents(contents) {
        var element = document.getElementsByClassName('btn btn-primary btn-sm');
        element.innerText = contents;
        initAceEditor()
      }
      
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
        });}