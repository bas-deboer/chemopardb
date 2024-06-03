/*global showAlert*/

function checkbox_selection () {
    $("#alignment_table > tbody > tr").click(function(event) {
        if (event.target.type !== "checkbox") {
            $(":checkbox", this).trigger("click");
            $(this).eq(0).toggleClass("alt_selected");
        }
        $(this).eq(0).toggleClass("alt_selected");
    });

    $("#select-all").click(function() {
        $(":checkbox", this).trigger("click");
        if ($(this).prop("checked")===true) {
            $(".alt:visible").prop("checked", true);
            $(".alt:visible").parent().parent().addClass("alt_selected");
        }
        if ($(this).prop("checked")===false) {
            $(".alt:visible").prop("checked", false);
            $(".alt:visible").parent().parent().removeClass("alt_selected");
        }
    });
}

function alignment_download () {
    $("#alignment_download_button").click(function() {
        var dataCSV = '';
        if ($(".alt_selected").length===0 || ($(".alt_selected").length===1 & $("#select-all").parent().parent().hasClass("alt_selected"))) {
            showAlert("No entries selected", "danger");
        }
        else {
            $("#alignment_table tr").each(function() {
                if ($(this).is(':first-child') || $(this).hasClass("alt_selected")) {
                    $(this).find("th").each(function(key, cell) {
                        dataCSV += $(cell).text()+",".repeat($(cell).prop("colSpan"));
                    })
                    $(this).find("td").each(function(key, cell) {
                        dataCSV += $(cell).text()+",";
                    })
                    dataCSV += "\n";
                }
            })
            dataCSV += consensus_symbols("#alignment_table");
            var hiddenElement = document.createElement('a');  
            hiddenElement.href = 'data:text/csv;charset=utf-8,' + encodeURI(dataCSV);  
            hiddenElement.target = '_blank';  
            hiddenElement.download = 'SH2DB_alignment.csv';  
            hiddenElement.click();
        }
    });
}

function mutation_symbols (table_id, domain) {
    var array = transposed_table(table_id, domain)
    var mutation_symbols = "<tr>";
    $(array).each(function(nums, vals) {
        if (vals.length<2 || (vals.includes("-") && vals.length==2) ) {
            mutation_symbols+= "<td></td>";
        }
        else {
             mutation_symbols+= "<td class='mutation'>*</td>";
        }
    })
    mutation_symbols+="</tr>"
    $(table_id + " tr." + domain + ":last").after(mutation_symbols);
}

function transposed_table (table_id, target_class) {
    var dataArray = [];
    $(table_id + " tr").each(function() {
        if ($(this).hasClass(target_class)) {
            var dataRow = [];
            $(this).find("td").each(function(key, cell) {
                dataRow.push($(cell).text());
            })
            dataArray.push(dataRow);
        }
    })
    if (dataArray.length===0) {
        return []
    }
    var transposedArray = [];
    for (var i=0; i<dataArray[0].length; i++) {
        var row = [];
        for (var j=0; j<dataArray.length; j++) {
            if (dataArray[j][i].length==1 && !row.includes(dataArray[j][i])) {
                row.push(dataArray[j][i]);
            }
        }
        transposedArray.push(row);
    }
    return transposedArray
}

function consensus_symbols (table_id) {
    var transposedArray = transposed_table (table_id, "alt_selected")
    var consensus_symbols = "";
    
    $(transposedArray).each(function(nums, vals) {
        skiprow = false;
        gap = false;
        var symbol = "";
        if (vals.length==1 && vals[0]!="-") {
            consensus_symbols += "*,";
        }
        else {
            if (vals.includes("-") || vals.length==0) {
                consensus_symbols += ",";
            }
            else {
                consensus_symbols += find_consensus_set(vals) + ",";
            }
        }
    })
    return consensus_symbols;
}

function find_consensus_set (array) {
    var colon_sets = [["S", "T", "A"],
                      ["N", "E", "Q", "K"],
                      ["N", "H", "Q", "K"],
                      ["N", "D", "E", "Q"],
                      ["Q", "H", "R", "K"],
                      ["M", "I", "L", "V"],
                      ["M", "I", "L", "F"],
                      ["H", "Y"],
                      ["F", "Y", "W"]];
    var dot_sets = [["C", "S", "A"],
                    ["A", "T", "V"],
                    ["S", "A", "G"],
                    ["S", "T", "N", "K"],
                    ["S", "T", "P", "A"],
                    ["S", "G", "N", "D"],
                    ["S", "N", "D", "E", "Q", "K"],
                    ["N", "D", "E", "Q", "H", "K"],
                    ["N", "E", "Q", "H", "R", "K"],
                    ["F", "V", "L", "I", "M"],
                    ["H", "F", "Y"]];
    
    for (j=0; j<colon_sets.length; j++) {
        var present = false;
        for (i=0; i<array.length; i++) {
            if (!colon_sets[j].includes(array[i])) {
                present = false;
                break;
            }
            else {
                present = true;
            }
        }
        if (present) {
            return ":";
        }
    }
    for (j=0; j<dot_sets.length; j++) {
        var present = false;
        for (i=0; i<array.length; i++) {
            if (!dot_sets[j].includes(array[i])) {
                present = false;
                break;
            }
            else {
                present = true;
            }
        }
        if (present) {
            return ".";
        }
    }
    return "";
}

function structure_download () {
    $("#structure_download_button").click(function() {
        var structures = [];
        if ($(".structure.alt_selected").length===0) {
            showAlert("No structure entries selected", "danger");
        }
        else {
            $(".structure.alt_selected").each(function() {
                if ($(this).hasClass('alphafold')) {
                    structures.push($(":nth-child(2)", this).text()+'_'+$(":nth-child(3)", this).text().slice(-1));
                }
                else {
                    structures.push($(":nth-child(3)", this).text());
                }
            });
            // StructureDownload(structures.join(','));
            window.location.href = '/structure/download?ids='+structures.join(",");
        }
    });
}

/* DEPRECATED, EVERY PAGE USES pymol_download_from_search NOW */
function pymol_download () {
    $("#pymol_download_button").click(function() {
        var structures = [];
        var residues = [];
        if ($(".structure.alt_selected").length===0) {
            showAlert("No structure entries selected", "danger");
        }
        else if ($(".residue_checkbox :checkbox:checked").length===0) {
            showAlert("No residues selected", "danger");
        }

        else {
            $(".structure.alt_selected").each(function() {
                structures.push($(":nth-child(3)", this).text());
            });
            $(".residue_checkbox :checkbox:checked").each(function() {
                residues.push($(this).attr('id'));
            });
            // PymolDownload(structures.join(','));
            window.location.href = '/structure/pymoldownload?ids='+structures.join(",")+'&residues='+residues.join(",");
        }
    });
}

function pymol_download_from_search () {
    $("#pymol_download_button").click(function() {
        var structures = [];
        var residues = [];
        var column_indices = [];
        if ($(".structure :checkbox:checked").length===0) {
            showAlert("No structure entries selected", "danger");
        }

        /*
        else if ($(".residue_checkbox :checkbox:checked").length===0) {
            showAlert("No residues selected", "danger");
        }
        */

        else {
            $(".residue_checkbox :checkbox:checked").each(function() {
                column_indices.push($(this).parent().parent().children().index($(this).parent()) );
            });

            $(".structure :checkbox:checked").each(function() {
                if($(this).parent().parent().hasClass('alphafold')) {
                    structures.push($(this).parent().parent().children().eq(1).text()+'_'+$(this).parent().parent().children().eq(2).text().slice(-1));
                }
                else {
                    structures.push($(this).parent().parent().children().eq(2).text() );
                }

                $.each(column_indices, function(key,value) {
                    if( $(this).parent().parent().children().eq(value).attr('title') ) {
                        residues.push($(this).parent().parent().children().eq(value).attr('title').split('\n')[0].replace(/[^0-9]/gi, '') );
                    }

                    else {
                        residues.push( '' );
                    }
                }.bind(this));

            });
            
            // PymolDownload(structures.join(','));
            window.location.href = '/structure/pymoldownload?ids='+structures.join(",")+'&residues='+residues.join(",");
        }
    });
}

function StructureDownload (structures) {
    $.ajax({
        'url': '/structure/download',
        'data': {"structures": structures},
        'type': 'GET',
        'success': function(response) {},
        'error': function(response) {}
    });
}

function table_filter () {
    fill_filter(".protein", "#protein_filter");
    fill_filter(".domain", "#domain_filter");
    for (var i=1; i<$(".data-row:first > .residue").length+1; i++) {
        fill_filter(".res_"+i.toString(), "#filter_"+i.toString())
    }
    $("select").change(function() {
        run_filter();
    });

}

function remove_filter_options () {
    $("#protein_filter").find("option").remove();
    $("#domain_filter").find("option").remove();
    for (var i=1; i<$(".data-row:first > .residue").length+1; i++) {
        $("#filter_"+i.toString()).find("option").remove();
    }
}

function update_filter_options () {
    remove_filter_options();
    fill_filter(".protein", "#protein_filter", true);
    fill_filter(".domain", "#domain_filter", true);
    for (var i=1; i<$(".data-row:first > .residue").length+1; i++) {
        fill_filter(".res_"+i.toString(), "#filter_"+i.toString(), true)
    }
}

function run_filter () {
    var selected_filters = {};
    $(".form-select").each(function (key, val) {
        var selected_value = $(val).find(":selected").text();
        if ($(val).find(":selected").text()!=="") {
            selected_filters[$(val).attr("id")] = selected_value;
        }
    })
    $(".data-row").addClass("hidden");
    var rows_to_show = [];
    $(".data-row").each(function (i, j) {
        if ($("#structure_toggle_button").hasClass("left_position") && $(j).hasClass("pdb")) {
            return;
        }
        var show = true;
        for (var [key, val] of Object.entries(selected_filters)) {
            if (key.startsWith("filter")) {
                key = key.replace("filter","res")
            }
            else {
                key = key.split('_')[0];
            }
            if ($(j).find("."+key).text()!==val) {
                show = false;
                return;
            }
        }
        if (show) {
            $(j).removeClass("hidden");
        }
    })
    table_numbers();
}

function fill_filter (td_class, filter_id, only_visible=false) {
    var options = [];
    $(td_class).each(function (key, val) {
        if (only_visible) {
            if (!$(val).parent().hasClass('hidden') && !options.includes(val.innerText)) {
                options.push(val.innerText);
            }
        }
        else {
            if (!options.includes(val.innerText)) {
                options.push(val.innerText);
            }
        }
    });
    $(filter_id).append(`<option value=""></option>`)
    $(options).each(function (key, val) {
        $(filter_id).append(`<option value="${val}">${val}</option>`);
    });
}

function sheinerman_button () {
    $("#sheinerman_button").click(function() {
        if (!$(this).hasClass("active")) {
            $(this).addClass("active");
            $(this).text("Unselect Sheinerman residues");
            var checkbox_status = true;
        }
        else {
            $(this).removeClass("active");
            $(this).text("Select Sheinerman residues");
            var checkbox_status = false;
        }
        
        $('.sheinerman').each(function (key, val) {
            $(".residue_checkbox").find("input").eq($(val).index()-3).prop("checked", checkbox_status)
        })
    })
}

function table_numbers () {
    var total = $('#alignment_table').find('tbody tr').length;
    var total_wt = $('#alignment_table').find('.alphafold').length;
    var visible = total-$('#alignment_table').find('tr.hidden').length;
    var visible_wt = total_wt-$('#alignment_table').find('.alphafold.hidden').length;
    var selected = $('#alignment_table').find('.alt_selected').length;
    var selected_wt = $('#alignment_table').find('.alt_selected.alphafold').length;
    $('#num_total').text(total);
    $('#num_total_wt').text(total_wt);
    $('#num_visible').text(visible);
    $('#num_visible_wt').text(visible_wt);
    $('#num_selected').text(selected);
    $('#num_selected_wt').text(selected_wt);
    
}
