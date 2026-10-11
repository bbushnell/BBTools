package prot;

import java.io.File;
import java.util.HashMap;

import fileIO.ByteFile;
import parse.LineParser1;
import structures.ByteBuilder;
import map.LongHashSet;

/**
 * Pinned sparse occurrence selection for original-corpus replay. Parent provenance
 * references are authenticated transitively by the independently supplied selector
 * pin; the extractor must verify the large parent ledgers once before producing it.
 * Replay verifies this file, the origin manifest, model, IDs, splits and ordinals.
 * @author Brian, Nilou
 */
final class MagQCReplaySelection {

	/** Loads a fresh, single-use selector before a replay spool or output is opened. */
	static MagQCReplaySelection load(String path, String pin, String model, String manifestPin){
		DigestSuffix.requireSuffix(pin, "replayselectsha80");
		if(!pin.equals(DigestSuffix.file(path))){throw new IllegalArgumentException("Replay selection hash mismatch: "+path);}
		final MagQCReplaySelection selection=new MagQCReplaySelection(new File(path).getAbsolutePath(), pin, model);
		final ByteFile input=ByteFile.makeByteFile(path, false);
		final LineParser1 parser=new LineParser1((byte)'\t');
		boolean schema=false, origin=false, named=false, columns=false;
		try{
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				parser.set(line);
				if(line.length==0){throw new IllegalArgumentException("Blank replay selection line");}
				if(line[0]=='#'){
					if(columns){throw new IllegalArgumentException("Late replay selection metadata");}
					if(parser.termEquals("#provenance", 0)){
						if(parser.terms()!=4){throw new IllegalArgumentException("Expected provenance kind/path/sha80");}
						final int kind=kind(parser.parseString(1));
						final String parent=parser.parseString(2), parentPin=parser.parseString(3);
						if(selection.parents[kind]!=null || !new File(parent).isAbsolute()){
							throw new IllegalArgumentException("Duplicate provenance or nonabsolute parent path");
						}
						DigestSuffix.requireSuffix(parentPin, "selection provenance");
						if((parser.termEquals("original_train", 1) || parser.termEquals("original_val", 1)) && !parentPin.equals(manifestPin)){
							throw new IllegalArgumentException("Original provenance must bind this selector's source manifest");
						}
						selection.parents[kind]=parent; selection.parentPins[kind]=parentPin;
						selection.parentKinds[kind]=parser.parseString(1);
					}else{
						if(parser.terms()!=2){throw new IllegalArgumentException("Expected two-field selection metadata");}
						if(parser.termEquals("#schema_version", 0) && !schema && parser.termEquals("magqc_replay_selection_v1", 1)){schema=true;}
						else if(parser.termEquals("#bin_manifest_sha80", 0) && !origin && parser.termEquals(manifestPin, 1)){origin=true;}
						else if(parser.termEquals("#bin_model", 0) && !named && parser.termEquals(model, 1)){named=true;}
						else{throw new IllegalArgumentException("Unknown, duplicate or mismatched replay selection metadata");}
					}
					continue;
				}
				if(!columns){
					if(!schema || !origin || !named || !java.util.Arrays.equals(line, COLUMNS)){
						throw new IllegalArgumentException("Missing selection bindings or invalid columns");
					}
					columns=true; continue;
				}
				if(parser.terms()!=5){throw new IllegalArgumentException("Expected five selection fields");}
				final String id=parser.parseString(0), split=parser.parseString(1);
				final long ordinal=parser.parseLong(2), parentRow=parser.parseLong(4);
				final int kind=kind(parser.parseString(3));
				if(id.isEmpty() || !split.equals(kind==0 ? "train" : "val") || ordinal<0 || ordinal>Integer.MAX_VALUE
						|| parentRow<0 || selection.parents[kind]==null
						|| !parser.termEquals(selection.parentKinds[kind], 3) || !selection.parentRows[kind].add(parentRow)){
					throw new IllegalArgumentException("Invalid or duplicate selection provenance for "+id);
				}
				if(selection.rows.put(id, new Row(split, (int)ordinal))!=null){
					throw new IllegalArgumentException("Duplicate selected bin ID: "+id);
				}
			}
		}finally{if(input.close()){throw new IllegalArgumentException("I/O failure reading replay selector");}}
		if(!columns || selection.rows.isEmpty()){throw new IllegalArgumentException("Empty replay selector");}
		if(!pin.equals(DigestSuffix.file(path))){throw new IllegalArgumentException("Replay selector changed while reading");}
		return selection;
	}

	private MagQCReplaySelection(String path_, String pin_, String model_){path=path_; pin=pin_; model=model_;}

	/** Validates a requested occurrence against its immutable original membership. */
	boolean retain(MagQCBinManifest.Bin bin, int ordinal){
		assert(bin!=null) : "Replay source reader must supply a parsed manifest row";
		return retain(bin.id,bin.split,ordinal);
	}

	/** Matches authenticated existing vector/BP identity without rebuilding source-bin geometry. */
	boolean retain(final String id, final String split, final int ordinal){
		final Row row=rows.get(id);
		if(row==null){return false;}
		if(row.seen || ordinal!=row.ordinal || !split.equals(row.split)){
			throw new IllegalArgumentException("Selected ID has wrong split/ordinal or appears twice: "+id);
		}
		row.seen=true; found++;
		return true;
	}

	/** Fails before inference if any requested ID was absent from the authenticated source. */
	void finish(){
		if(found!=rows.size()){
			for(java.util.Map.Entry<String,Row> entry:rows.entrySet()){
				if(!entry.getValue().seen){throw new IllegalArgumentException("Selected ID missing from manifest: "+entry.getKey());}
			}
			throw new IllegalStateException("Selection accounting mismatch");
		}
	}

	/** Pins output to the exact selector, including the parent-row mapping it carries. */
	void appendHeader(ByteBuilder out){
		assert(found==rows.size()) : "No output may advertise an incompletely matched selector";
		out.append("#replay_selection_file\t").append(path).nl();
		out.append("#replay_selection_sha80\t").append(pin).nl();
		out.append("#replay_source_model\t").append(model).nl();
	}

	private static int kind(String name){
		if("old5m".equals(name) || "original_train".equals(name)){return 0;}
		if("original400k".equals(name) || "original_val".equals(name)){return 1;}
		throw new IllegalArgumentException("Unknown selection provenance kind: "+name);
	}

	private static final class Row {
		Row(String split_, int ordinal_){split=split_; ordinal=ordinal_;}
		final String split;
		final int ordinal;
		boolean seen;
	}
	private final String path, pin, model;
	private final HashMap<String,Row> rows=new HashMap<String,Row>();
	private final String[] parents=new String[2], parentPins=new String[2], parentKinds=new String[2];
	private final LongHashSet[] parentRows={new LongHashSet(), new LongHashSet()};
	private int found;
	private static final byte[] COLUMNS="bin_id\tsplit\tsource_ordinal\tprovenance\tprovenance_row".getBytes(java.nio.charset.StandardCharsets.US_ASCII);
}
