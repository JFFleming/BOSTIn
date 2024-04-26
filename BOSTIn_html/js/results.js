var resultsPage;
function resetForm() {
    document.getElementById('faFile').value = "";
}

function bostin() {
    console.log("BOSTIn function called.");
    var input = document.getElementById("faFile");
    var files = input.files;
    if (files.length == 0) {
        alert("Please upload an aligned fasta file!");
        return;
    }
    var inputFile = input.value;

    // Get the selected values
    var selectElement = document.getElementById('type');
    var dtype = selectElement.value;
    if (dtype === "none") {
        alert("Please select a datatype for your alignment - either nucleotide or amino acids");
        return;
    }

    var blhCheckbox = document.getElementById('BLH').checked;
    var siteSatCheckbox = document.getElementById('SiteSat').checked;
    var compHetCheckbox = document.getElementById('CompHet').checked;

    if (!blhCheckbox && !siteSatCheckbox && !compHetCheckbox) {
        alert("No assessments have been checked. Please select at least one option.");
        return;
    }
    if (siteSatCheckbox && dtype.match("nucl")) {
        alert("Sorry, Site Saturation for DNA is not yet implemented")
	    return;
	}

    var resultsPage = window.open("results.html");

    console.log("Opening results page...");
    // Perform actions based on user selections
	resultsPage.onload = function() {
		console.log("function working");
		var resultsTable = resultsPage.document.getElementById("resultsTable");

	    if (blhCheckbox) {
	    	console.log("blh checkbox");
	        // Perform action for BLH checkbox checked
	        // For example: calculate BLH score
	        var blhResult = "boop";
	        addResultToTable("Branch Length Heterogeneity (LB-Score)", blhResult, resultsPage);
	    }
	    if (siteSatCheckbox && dtype.match("prot")) {
	    	console.log("sat checkbox");
            // Perform action for SiteSat checkbox checked
            // For example: calculate SiteSat score
            var siteSatResult = "foo";
            addResultToTable("Site Saturation (DE-Score)", siteSatResult, resultsPage);
        }
	    if (compHetCheckbox) {
	    	console.log("ch checkbox");
	    	// Perform action for CompHet checkbox checked
	    	// For example: calculate CompHet score
	    	var compHetResult = "bar";
	    	addResultToTable("Compositional Heterogeneity (nRCFV)", compHetResult, resultsPage);
	    }
    }
}
    
function addResultToTable(resultType, resultValue) {
    // Access the table in the results page
    var tableBody = resultsPage.document.querySelector('#resultsTable tbody');

    // Create a new row
    var newRow = resultsPage.document.createElement('tr');

    // Create cells for result type and result value
    var typeCell = resultsPage.document.createElement('td');
    typeCell.textContent = resultType;
    var valueCell = resultsPage.document.createElement('td');
    valueCell.textContent = resultValue;

    // Append cells to the row
    newRow.appendChild(typeCell);
    newRow.appendChild(valueCell);

    // Append the row to the table
    tableBody.appendChild(newRow);
}