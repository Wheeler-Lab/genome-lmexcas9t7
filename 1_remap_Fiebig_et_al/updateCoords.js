var fs=require("fs");
var path=require("path");

json=require(path.resolve(process.cwd(), process.argv[2]));

function mostFrequent(children, type) {
	var samples={"Promastigote": {}, "Intracellular_Amastigote": {}, "Axenic_Amastigote": {}};
	var maxFreq=0;
	var index=-1;
	for (child in children) {
		if (children[child].type==type) {
			for (sample in samples) {
				if (children[child].params[sample]) {
					if (parseInt(children[child].params[sample])>maxFreq) {
						maxFreq=children[child].params[sample];
						index=child;
					}
				}
			}
		}
	}
	return index;
}

function getTypeIndex(children, type) {
	for (var child in children) if (children[child].type==type) return child;
	return -1;
}

function getTypeExtent(children, type) {
	var res={};
	for (child in children) if (children[child].type==type) {
		if (!res.start) res.start=children[child].start;
		if (!res.end) res.end=children[child].end;
		res.start=Math.min(children[child].start, res.start);
		res.end=Math.max(children[child].end, res.end);
	}
	return res;
}

//Update gene IDs to standardised naming
//<geneID> unchanged
//mRNA named <geneID>.X
//exons named exon_<geneID>-EX
//CDSs named <geneID>-pX-CDSX
//Proteins named <geneID>-pX
for (var i in json) {
	var dir=json[i].dir;
	if (json[i].type=="gene") {
		if (json[i].children) {
			var mrnaindex=1;
			if (dir=="+") {
				var j=0;
			} else if (dir=="-") {
				var j=json[i].children.length-1;
			}
			while (j>=0 && j<json[i].children.length) {
				if (json[i].children[j].type=="mRNA") {
					var mrnaid=json[i].id+"."+mrnaindex
					json[i].children[j].id=mrnaid;
					json[i].children[j].params.ID=mrnaid;
					json[i].children[j].params.Parent=json[i].id;
					mrnaindex++;
					if (json[i].children[j].children) {
						var exonindex=1;
						var cdsindex=1;
						if (dir=="+") {
							var k=0;
						} else if (dir=="-") {
							var k=json[i].children[j].children.length-1;
						}
						while (k>=0 && k<json[i].children[j].children.length) {
							if (json[i].children[j].children[k].type=="exon") {
								var exonid="exon_"+json[i].id+"-E"+exonindex;
								json[i].children[j].children[k].id=exonid;
								json[i].children[j].children[k].params.ID=exonid;
								json[i].children[j].children[k].params.Parent=json[i].children[j].id;
								exonindex++;
							} else if (json[i].children[j].children[k].type=="CDS") {
								var cdsid=json[i].id+"-p"+cdsindex+"-CDS"+cdsindex;
								json[i].children[j].children[k].id=cdsid;
								json[i].children[j].children[k].params.ID=cdsid;
								json[i].children[j].children[k].params.Parent=json[i].children[j].id;
								var protid=json[i].id+"-p"+cdsindex;
								json[i].children[j].children[k].params.protein_source_id=protid;
								cdsindex++;
							}
							if (json[i].dir=="+") {
								k++;
							} else if (dir=="-") {
								k--;
							}
						}
					}
				}
				if (json[i].dir=="+") {
					j++;
				} else if (dir=="-") {
					j--;
				}
			}
		}
	}
}

//Add UTRs and correct gene coordinates
//UTRs named utr_<geneID>-UX
for (var entry in json) {
	if (json[entry].type=="gene") {
		for (var child=0; child<json[entry].children.length; child++) {
			if (json[entry].children[child].contig!=json[entry].contig) {
				//Remove children which have not ended up on the correct contig
				json[entry].children.splice(child, 1);
				child--;
			}
		}
		for (var child=0; child<json[entry].children.length; child++) {
			//Determine gene extent using most frequent PAS and SLAS
			var pasIndex=mostFrequent(json[entry].children, "PAS");
			var mrna=getTypeIndex(json[entry].children, "mRNA");
			if (mrna!=-1) {
				var cdsrange=getTypeExtent(json[entry].children[mrna].children, "CDS");
				if (pasIndex!=-1) {
					if (json[entry].dir=="+") {
						// Check if PAS is inside the CDS, if yes, skip defining a 5' UTR and leave gene start alone.
						if (json[entry].children[pasIndex].start > cdsrange.end) {
							json[entry].end=json[entry].children[pasIndex].end;
							// Update mRNA end to the end of the 3'UTR
							json[entry].children[mrna].end = json[entry].end;
							json[entry].children[mrna].children.push({
								type: "three_prime_UTR",
								id: "utr_"+json[entry].children[mrna].id+"-U1",
								contig: json[entry].contig,
								start: cdsrange.end+1,
								end: json[entry].end,
								dir: json[entry].dir,
								params: {
									ID: "utr_"+json[entry].children[mrna].id+"-U1",
									Parent: json[entry].children[mrna].id
								}
							});
						}
					} else if (json[entry].dir=="-") {
						// Check if PAS is inside the CDS, if yes, skip defining a 5' UTR and leave gene start alone.
						if (json[entry].children[pasIndex].start < cdsrange.start) {
							json[entry].start=json[entry].children[pasIndex].start;
							// Update mRNA end to the end of the 3'UTR
							json[entry].children[mrna].start = json[entry].start;
							json[entry].children[mrna].children.push({
								type: "three_prime_UTR",
								id: "utr_"+json[entry].children[mrna].id+"-U11",
								contig: json[entry].contig,
								start: json[entry].start,
								end: cdsrange.start-1,
								dir: json[entry].dir,
								params: {
									ID: "utr_"+json[entry].children[mrna].id+"-U1",
									Parent: json[entry].children[mrna].id
								}
							});
						}
					}
				}
				var slasIndex=mostFrequent(json[entry].children, "SLAS");
				if (slasIndex!=-1) {
					if (json[entry].dir=="+") {
						// Check if SLAS is inside the CDS, if yes, skip defining a 5' UTR and leave gene start alone.
						if (json[entry].children[slasIndex].start < cdsrange.start) {
							json[entry].start=json[entry].children[slasIndex].start;
							// Update mRNA start to the start of the 5'UTR
							json[entry].children[mrna].start = json[entry].start;
							json[entry].children[mrna].children.push({
								type: "five_prime_UTR",
								id: "utr_"+json[entry].children[mrna].id+"-U2",
								contig: json[entry].contig,
								start: json[entry].start,
								end: cdsrange.start-1,
								dir: json[entry].dir,
								params: {
									ID: "utr_"+json[entry].children[mrna].id+"-U2",
									Parent: json[entry].children[mrna].id
								}
							});
						}
					} else if (json[entry].dir=="-") {
						// Check if SLAS is inside the CDS, if yes, skip defining a 5' UTR and leave gene start alone.
						if (json[entry].children[slasIndex].end > cdsrange.end) {
							json[entry].end=json[entry].children[slasIndex].end;
							// Update mRNA start to the start of the 5'UTR
							json[entry].children[mrna].end = json[entry].end;
							json[entry].children[mrna].children.push({
								type: "five_prime_UTR",
								id: "utr_"+json[entry].children[mrna].id+"-U2",
								contig: json[entry].contig,
								start: cdsrange.end+1,
								end: json[entry].end,
								dir: json[entry].dir,
								params: {
									ID: "utr_"+json[entry].children[mrna].id+"-U2",
									Parent: json[entry].children[mrna].id
								}
							});
						}
					}
				}
			}
		}
	}
}

function recursiveRemoveParams(json) {
	var validParams={"ID": {}, "Parent": {}, "protein_source_id": {}, "Promastigote": {}, "Intracellular_Amastigote": {}, "Axenic_Amastigote": {}, "cdsHistory": {}};
	for (var i in json) {
		if (json[i].children) json[i].children=recursiveRemoveParams(json[i].children);
		if (json[i].params) for (param in json[i].params) if (!validParams[param]) delete json[i].params[param];
	}
	return json;
}

//Recursively clear unused parameters
json=recursiveRemoveParams(json);

console.log(JSON.stringify(json, "", "\t"));
