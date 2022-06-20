// Enable popover
var popoverTriggerList = [].slice.call(document.querySelectorAll('[data-bs-toggle="popover"]'))
var popoverList = popoverTriggerList.map(function (popoverTriggerEl) {
  return new bootstrap.Popover(popoverTriggerEl)
})
// Enable popover end

//tf

const expVal = document.getElementById("expression");
const trshRadio = document.querySelectorAll(".trsh-radio");

expVal.addEventListener("click", function(){
  trshRadio.disabled  = true;
})

//tf end