document.addEventListener("DOMContentLoaded", function () {
    const selectedDescriptors = new Set();

    window.showTab = function showTab(event, tabId) {
        var contents = document.querySelectorAll(".tab-content");
        contents.forEach(content => content.style.display = "none");

        var tabs = document.querySelectorAll(".tab");
        tabs.forEach(tab => tab.classList.remove("active"));

        var activeContent = document.getElementById(tabId);
        if (activeContent) {
            activeContent.style.display = "block";
        }
        event.target.classList.add("active");

        if (tabId !== "search-tab") {
            var searchCheckboxes = document.querySelectorAll('#searchResults input[type="checkbox"]');
            searchCheckboxes.forEach(function (sCheckbox) {
                var mainCheckbox = document.querySelector(
                    '.tab-content:not(#search-tab) .checkbox-container input[value="' + sCheckbox.value + '"]'
                );
                if (mainCheckbox) {
                    sCheckbox.checked = mainCheckbox.checked;
                }
            });
        }
    };

    function toggleSelectDeselect() {
        var selectAllButton = document.getElementById("selectAllButton");
        var checkboxes = document.querySelectorAll('#common-descriptors .checkbox-container input[type="checkbox"]');

        if (selectAllButton.innerHTML === "Select All") {
            checkboxes.forEach(cb => {
                cb.checked = true;
                selectedDescriptors.add(cb.value);
            });
            selectAllButton.innerHTML = "Deselect All";
        } else {
            checkboxes.forEach(cb => {
                cb.checked = false;
                selectedDescriptors.delete(cb.value);
            });
            selectAllButton.innerHTML = "Select All";
        }

        updateSearchResultCheckboxes();
    }

    function updateSearchResultCheckboxes() {
        var searchCheckboxes = document.querySelectorAll('#searchResults input[type="checkbox"]');
        searchCheckboxes.forEach(function (sCheckbox) {
            if (selectedDescriptors.has(sCheckbox.value)) {
                sCheckbox.checked = true;
            } else {
                sCheckbox.checked = false;
            }
        });
    }

    function filterDescriptors() {
        var input = document.getElementById("searchBar");
        var filter = input.value.toLowerCase().trim();
        var searchResults = document.getElementById("searchResults");
        searchResults.innerHTML = "";  // Clear previous search results

        if (!filter) return;  // If no input, return

        var addedDescriptors = new Set();

        for (var descriptor in descriptorsAndSynonyms) {
            var synonyms = descriptorsAndSynonyms[descriptor];
            var allTerms = [descriptor].concat(synonyms);

            for (var i = 0; i < allTerms.length; i++) {
                var term = allTerms[i].toLowerCase().trim();
                if (term.includes(filter)) {
                    if (addedDescriptors.has(descriptor)) {
                        continue;
                    }

                    addedDescriptors.add(descriptor);

                    var label = document.createElement("label");
                    var inputCheckbox = document.createElement("input");
                    inputCheckbox.type = "checkbox";
                    inputCheckbox.name = "displayOptions";
                    inputCheckbox.value = descriptor;
                    
                    if (selectedDescriptors.has(descriptor)) {
                        inputCheckbox.checked = true;
                    }

                    label.appendChild(inputCheckbox);

                    var displayTerm = allTerms[i];
                    if (displayTerm === descriptor) {
                        label.appendChild(document.createTextNode(" " + descriptor));
                    } else {
                        label.appendChild(document.createTextNode(" " + descriptor + " (" + displayTerm + ")"));
                    }

                    inputCheckbox.addEventListener("change", function (e) {
                        var value = this.value;
                        if (this.checked) {
                            selectedDescriptors.add(value);
                        } else {
                            selectedDescriptors.delete(value);
                        }

                        var mainCheckbox = document.querySelector(
                            '.tab-content:not(#search-tab) .checkbox-container input[value="' + value + '"]'
                        );
                        if (mainCheckbox) {
                            mainCheckbox.checked = this.checked;
                        }
                    });

                    searchResults.appendChild(label);
                    break;
                }
            }
        }

        updateSearchResultCheckboxes();
    }

    window.jsmeOnLoad = function () {
        if (!window.jsmeApplet) {
            window.jsmeApplet = new JSApplet.JSME("jsme_container", "380px", "340px");
            window.jsmeApplet.setAfterStructureModifiedCallback(function () {
                var smiles = window.jsmeApplet.smiles();
                var formattedSmiles = smiles.replace(/\n/g, ",");
                document.getElementById("smiles").value = formattedSmiles;
            });
        }
    };

    // Open or close the molecule drawer
    window.openMoleculeDrawer = function () {
        var jsmeContainer = document.getElementById("jsme_container");
        var drawButton = document.querySelector(".draw-button");

        if (jsmeContainer.style.display === "none") {
            jsmeContainer.style.display = "block";
            jsmeOnLoad();  // Load JSME if not already loaded
            drawButton.classList.add("draw-button-active");
            drawButton.textContent = "Close";
        } else {
            jsmeContainer.style.display = "none";
            drawButton.classList.remove("draw-button-active");
            drawButton.textContent = "Draw Molecules";
        }
    };

    var searchBar = document.getElementById("searchBar");
    if (searchBar) {
        let debounceTimeout;
        searchBar.addEventListener("keyup", function () {
            clearTimeout(debounceTimeout);
            debounceTimeout = setTimeout(filterDescriptors, 300);  // 300ms delay
        });
    }

    document.querySelectorAll('.tab').forEach(function (tab) {
        tab.addEventListener('click', function (event) {
            const tabId = tab.getAttribute('onclick').match(/'(.*)'/)[1];
            showTab(event, tabId);
        });
    });

    var csvFileUpload = document.getElementById("csvFileUpload");
    if (csvFileUpload) {
        csvFileUpload.addEventListener("change", function () {
            if (this.value) {
                this.form.submit();
            }
        });
    }

    var selectAllButton = document.getElementById("selectAllButton");
    if (selectAllButton) {
        selectAllButton.addEventListener("click", toggleSelectDeselect);
    }

});
