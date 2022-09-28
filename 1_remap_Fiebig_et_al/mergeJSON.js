var fs=require("fs");
var path=require("path");

//Original JSON
json=require(path.resolve(process.cwd(), process.argv[2]));
//Additional JSON
addjson=require(path.resolve(process.cwd(), process.argv[3]));

function recursiveFlagCds(json) {
	for (var entry in json) {
		if (json[entry].type=="CDS") addjson[entry].params.cdsHistory="Added_FiebigEtAl2015";
		if (json[entry].children) json[entry].children=recursiveFlagCds(json[entry].children);
	}
	return json;
}

addjson=recursiveFlagCds(addjson);

function recursiveAddEntries(json, addjson) {
	for (var addentry in addjson) {
		var entryindex=-1;
		for (var entry in json) if (json[entry].id==addjson[addentry].id) entryindex=entry;
		if (entryindex==-1) {
			//If not already in the json
			json.push(addjson[addentry]);
		} else {
			//Otherwise, check and merge children
			if (json[entryindex].children && addjson[addentry].children) {
				//json[entryindex].children.push(addjson[addentry].children);
				json[entryindex].children=recursiveAddEntries(json[entryindex].children, addjson[addentry].children);
			}
		}
	}
	return json;
}

json=recursiveAddEntries(json, addjson);

//Remove entries where root level entries lack coordinates or a contig
for (var i=0; i<json.length; i++) {
	if (!json[i].start || !json[i].end || !json[i].contig || !json[i].dir) {
		console.error(json[i].id+" lacks coordinates, removed");
		json.splice(i, 1);
		i--;
	}
}

function recursiveSort(json) {
	var maxcoord=0;
	for (var entry in json) if (json[entry].end>maxcoord) maxcoord=json[entry].end;
	json.sort();
	json.sort(function(a, b) {
		if (!a.contig) console.error(a);
		if (!b.contig) console.error(b);
		var contigcompare=a.contig.toString().localeCompare(b.contig.toString());
		return a.start-b.start+contigcompare*maxcoord;
	});
}

recursiveSort(json);

console.log(JSON.stringify(json, "", "\t"));