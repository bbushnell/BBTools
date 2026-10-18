package prok;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;
import java.util.TreeMap;
import java.util.TreeSet;

import dna.AminoAcid;
import dna.GeneticCode;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import gff.GffLine;
import shared.Tools;
import structures.ByteBuilder;

/** Versioned reference concordance, not a biological true/false gene classifier.
 * The normative schema is docs/bbgap/STOP_EVALUATION.md.
 * @author Sayu
 */
public final class ProteinStopEvaluation {

	static void run(String[] args){
		String fasta=null, ref=null, baseline=null, candidate=null, out=null;
		int table=11;
		final HashSet<String> seen=new HashSet<String>();
		for(String arg : args){
			final int eq=arg.indexOf('=');
			require(eq>0 && eq<arg.length()-1, "Expected nonempty flag=value: "+arg);
			final String key=arg.substring(0, eq).toLowerCase(Locale.ROOT), value=arg.substring(eq+1);
			require(seen.add(key), "Duplicate option: "+key);
			if(key.equals("mode")){require(value.equalsIgnoreCase(SCHEMA), "Unknown evaluation mode: "+value);}
			else if(key.equals("fasta")){fasta=value;}
			else if(key.equals("ref")){ref=value;}
			else if(key.equals("baseline")){baseline=value;}
			else if(key.equals("candidate")){candidate=value;}
			else if(key.equals("out")){out=value;}
			else if(key.equals("table")){table=Integer.parseInt(value);}
			else{throw new IllegalArgumentException("Unknown stop evaluation option: "+key);}
		}
		require(fasta!=null && ref!=null && out!=null, "stops-v1 requires fasta= ref= out=PREFIX");
		require((baseline==null)==(candidate==null), "Specify both baseline and candidate, or neither for reference preflight");
		GeneticCode.forTable(table);
		final String[] outputs=baseline==null ? new String[]{out+".reference.tsv", out+".summary.tsv", out+".contigs.tsv"}
			: new String[]{out+".reference.tsv", out+".baseline.tsv", out+".candidate.tsv", out+".stops.tsv", out+".summary.tsv", out+".contigs.tsv"};
		require(Tools.testInputFiles(false, true, fasta, ref, baseline, candidate), "Unreadable stop evaluation input");
		require(Tools.testOutputFiles(false, false, false, outputs), "Outputs must be new writable files");
		final HashMap<String,byte[]> sequences=AnnotateGenes.readScaffolds(fasta);
		final Dataset reference=normalize(readGff(ref), sequences, table, true);
		final Dataset before=baseline==null ? null : normalize(readGff(baseline), sequences, table, false);
		final Dataset after=candidate==null ? null : normalize(readGff(candidate), sequences, table, false);
		if(before!=null){applyMasks(before, reference); applyMasks(after, reference);}
		// No output writer opens until all three inputs pass validation.
		writeContigs(out+".contigs.tsv", sequences, reference.circular);
		writeRows(out+".reference.tsv", reference);
		if(before!=null){
			writeRows(out+".baseline.tsv", before); writeRows(out+".candidate.tsv", after);
			writeStops(out+".stops.tsv", reference, before, after);
		}
		final ByteStreamWriter writer=new ByteStreamWriter(out+".summary.tsv", false, false, false); writer.start();
		try{
			writer.println("schema\tarm\tmetric\tvalue");
			writeSummary(writer, "reference", reference, reference);
			if(before!=null){
				writeSummary(writer, "baseline", before, reference); writeSummary(writer, "candidate", after, reference);
				writeDeltaSummary(writer, reference, before, after);
			}
		}finally{require(!writer.poisonAndWait(), "Failed writing stop summary");}
		System.err.println("STOP_EVALUATION_PASS schema="+SCHEMA+" reference_only="+(before==null));
	}

	/** Native line reader plus GffLine, with a strict nine-column guard before its legacy parser. */
	public static ArrayList<GffLine> readGff(String path){
		final ArrayList<GffLine> rows=new ArrayList<GffLine>();
		final ByteFile file=ByteFile.makeByteFile(path, true);
		final boolean oldSource=GffLine.paseSource;
		// The native parser otherwise discards column 2, defeating producer-specific adapters.
		GffLine.paseSource=true;
		try{
			for(byte[] line=file.nextLine(); line!=null; line=file.nextLine()){
				if(line.length==0){continue;}
				if(line[0]=='#'){
					require(!new String(line, java.nio.charset.StandardCharsets.UTF_8).equals("##FASTA"), "Embedded FASTA is unsupported; supply fasta= separately");
					continue;
				}
				int tabs=0; for(byte b : line){if(b=='\t'){tabs++;}}
				require(tabs==8, "Expected nine GFF columns in "+path);
				rows.add(new GffLine(line));
			}
		}finally{GffLine.paseSource=oldSource; require(!file.close(), "Failed reading GFF: "+path);}
		return rows;
	}

	public static Dataset normalize(ArrayList<GffLine> input, HashMap<String,byte[]> sequences, int fallback, boolean reference){
		final Dataset data=new Dataset(reference);
		if(reference){
			final HashMap<String,String> markers=new HashMap<String,String>();
			for(GffLine g : input){
				final String marker=attr(g, "Is_circular");
				if(!g.type.equals("region") || marker==null){continue;}
				require(marker.equals("true") || marker.equals("false"), "Invalid circular region marker: "+coordinate(g));
				final String old=markers.put(g.seqid, marker);
				require(old==null || old.equals(marker), "Conflicting circular region markers: "+g.seqid);
				if(marker.equals("true")){
					validateBounds(g, sequences);
					require(g.start==1 && g.stop==sequences.get(g.seqid).length, "Circular region must match full FASTA length: "+coordinate(g));
					data.circular.add(g.seqid);
				}
			}
		}
		final HashMap<String,GffLine> parents=new HashMap<String,GffLine>();
		final HashSet<String> ambiguousParents=new HashSet<String>();
		final HashMap<String,Group> groups=new HashMap<String,Group>();
		for(GffLine g : input){
			if(!g.type.equals("CDS")){
				final String id=attr(g, "ID");
				if(id!=null && parents.put(id, g)!=null){ambiguousParents.add(id);}
			}
		}
		int index=0;
		for(GffLine g : input){
			index++;
			if(!g.type.equals("CDS")){continue;}
			final boolean wraps=validateInterval(g, sequences, data.circular);
			final Row row=new Row(g, index); data.rows.add(row); data.rawCds++;
			row.circular=wraps;
			final String id=attr(g, "ID"), parent=attr(g, "Parent");
			// Identical protein accessions at different loci must not merge.
			final String group=id==null && parent==null ? "row:"+index : g.seqid+"\t"+g.strand+"\t"+parent+"\t"+id;
			Group members=groups.get(group);
			if(members==null){members=new Group(); groups.put(group, members);}
			members.rows.add(row); row.group=members;
			String reason=null;
			if(parent!=null){
				if(parent.indexOf(',')>=0){reason="multiple_parents";}
				else if(ambiguousParents.contains(parent)){reason="ambiguous_parent";}
				else{
					final GffLine p=parents.get(parent);
					if(p==null){reason="missing_parent";}
					else{
						// Validate before linkage exclusions; malformed virtual parents must never bypass bounds checks.
						final boolean parentWraps=validateInterval(p, sequences, data.circular);
						final int n=sequences.get(p.seqid).length;
						final boolean contains=(p.start<=g.start && p.stop>=g.stop)
							|| (parentWraps && g.stop<=p.stop-n);
						if(!p.seqid.equals(g.seqid) || (p.strand<=1 && p.strand!=g.strand) || !contains){reason="inconsistent_parent";}
						else{
							if(parentWraps){row.circularParent=p;}
							if(attr(p, "Parent")!=null){reason="parent_chain";}
							else{row.parent=p;}
						}
					}
				}
			}
			row.reason=reason;
		}
		for(Row row : data.rows){
			final GffLine g=row.g;
			final String declared=attr(g, "transl_table");
			row.table=declared==null ? fallback : Integer.parseInt(declared);
			final Group group=row.group;
			if(row.reason==null && row.circular){row.reason="circular_origin";}
			if(row.reason==null && row.circularParent!=null){row.reason="circular_parent";}
			// Equal-coordinate repeats are duplicates, not a multipart gene.
			final HashSet<String> coords=new HashSet<String>();
			for(Row member : group.rows){coords.add(coordinate(member.g));}
			if(row.reason==null && coords.size()>1){row.reason="multipart";}
			if(pseudo(g) || (row.parent!=null && pseudo(row.parent))){row.reason="pseudo";}
			if(row.reason==null){row.reason=classify(row, sequences.get(g.seqid), reference);}
			if(row.reason==null){row.reason="primary";}
			if(!reference && (row.reason.equals("missing_terminal_stop") || row.reason.equals("internal_stop") || row.reason.equals("incomplete_terminal_codon"))){
				throw new IllegalArgumentException("Invalid unmarked predicted CDS: "+coordinate(g)+" "+row.reason);
			}
		}
		if(reference){
			for(GffLine g : input){
				if(!g.type.equals("CDS") && pseudo(g)){
					validateInterval(g, sequences, data.circular); addMask(data, g, g.strand, sequences);
					final Row row=new Row(g, 0); row.reason="pseudo_mask"; data.maskRows.add(row);
				}
			}
			for(Group group : groups.values()){
				int low=Integer.MAX_VALUE, high=0; boolean excluded=false, circular=false;
				for(Row row : group.rows){
					low=Math.min(low, row.g.start); high=Math.max(high, row.g.stop); excluded|=!row.primary();
					circular|=row.circular || row.circularParent!=null;
				}
				if(excluded){
					if(circular){
						// A physical min/max hull across the origin would mask unrelated middle-of-contig genes.
						final HashSet<String> masked=new HashSet<String>();
						for(Row row : group.rows){
							final GffLine g=row.circularParent==null ? row.g : row.circularParent;
							// Parent supplies the arc, but a known child strand must not become a both-strand group mask.
							if(masked.add(coordinate(g))){addMask(data, g, row.g.strand, sequences);}
						}
					}else{final GffLine first=group.rows.get(0).g; data.masks.add(new Mask(first.seqid, first.strand, low, high));}
				}
			}
		}
		data.index();
		assert(data.rawCds==data.rows.size()) : "Every raw CDS must survive normalization with an explicit stratum";
		return data;
	}

	private static String classify(Row row, byte[] bases, boolean reference){
		final GffLine g=row.g;
		if(g.strand!=GffLine.PLUS && g.strand!=GffLine.MINUS){return "unknown_strand";}
		if(attr(g, "exception")!=null || attr(g, "transl_except")!=null || (row.parent!=null &&
			(attr(row.parent, "exception")!=null || attr(row.parent, "transl_except")!=null))){return "translation_exception";}
		// These input conventions need their own verified adapter, not speculative coordinate conversion.
		if(attr(g, "start_range")!=null || attr(g, "end_range")!=null || (row.parent!=null &&
			(attr(row.parent, "start_range")!=null || attr(row.parent, "end_range")!=null))){return "uncertain_range";}
		g.partial();
		final String partial=attr(g, "partial"), pp=row.parent==null ? null : attr(row.parent, "partial");
		if(row.parent!=null){row.parent.partial();}
		if(pp!=null && !pp.equals("false") && !pp.equals("00")){return "parent_partial";}
		if("true".equals(partial)){return "partial_side_unknown";}
		final boolean left=partial!=null && partial.length()==2 && partial.charAt(0)=='1';
		final boolean right=partial!=null && partial.length()==2 && partial.charAt(1)=='1';
		row.five=g.strand==GffLine.PLUS ? left : right;
		row.three=g.strand==GffLine.PLUS ? right : left;
		if(row.three){return "partial_3prime";}
		if(g.phase<0 || g.phase>2){return "unknown_phase";}
		if(row.table!=4 && row.table!=11 && row.table!=25){return "unsupported_table";}
		final GeneticCode code=GeneticCode.forTable(row.table);
		final int length=g.stop-g.start+1, translated=length-g.phase;
		if(translated<6){return "empty_cds";}
		if(translated%3!=0){return "incomplete_terminal_codon";}
		final int terminal=codon(g, bases, length-3);
		// Orf.appendGff emits stop:NNN for internal ambiguous boundaries, while
		// Orf.partialMask intentionally marks only physical contig-edge truncation.
		if(!reference && unknownNativeStop(g)){
			require(terminal<0, "BBTools stop:NNN contradicts an unambiguous FASTA triplet: "+coordinate(g));
			return "synthetic_stop";
		}
		if(!code.isStop(terminal)){return "missing_terminal_stop";}
		for(int offset=g.phase; offset<length-3; offset+=3){
			final int codon=codon(g, bases, offset);
			if(codon<0){return "ambiguous_codon";}
			if(code.isStop(codon)){return "internal_stop";}
		}
		return null;
	}

	private static boolean unknownNativeStop(GffLine g){
		final String a=g.attributes;
		if(!"BBTools".equals(g.source) || a==null || !a.startsWith("CDS,")){return false;}
		final int semicolon=a.indexOf(';'), limit=semicolon<0 ? a.length() : semicolon;
		boolean seen=false, unknown=false;
		for(int from=4; from<limit; ){
			int end=a.indexOf(',', from); if(end<0 || end>limit){end=limit;}
			if(end-from>=5 && a.regionMatches(from, "stop:", 0, 5)){
				require(!seen, "Duplicate native stop field: "+coordinate(g)); seen=true;
				unknown=end-from==8 && a.regionMatches(from, "stop:NNN", 0, 8);
			}
			from=end+1;
		}
		return unknown;
	}

	private static int codon(GffLine g, byte[] bases, int offset){
		assert(offset>=0 && offset+3<=g.stop-g.start+1) : "Endpoint validation must read a complete oriented codon";
		if(g.strand==GffLine.PLUS){return GeneticCode.codon(bases, g.start-1+offset);}
		int value=0;
		for(int i=0; i<3; i++){
			final int b=bases[g.stop-1-offset-i];
			final int n=b<0 || b>=AminoAcid.baseToNumber.length ? -1 : AminoAcid.baseToNumber[b];
			if(n<0 || b=='U' || b=='u'){return -1;}
			value=(value<<2)|(3-n);
		}
		return value;
	}

	public static void applyMasks(Dataset calls, Dataset ref){
		assert(!calls.reference && ref.reference) : "Only fixed reference exclusions may mask predictions";
		for(Row row : calls.rows){
			if(!row.primary() || ref.primary.containsKey(row.stop)){continue;}
			for(Mask mask : ref.masks){
				if(mask.contains(row.g)){row.reason="reference_exception_mask"; break;}
			}
		}
		calls.index();
	}

	private static void validateBounds(GffLine g, HashMap<String,byte[]> sequences){
		final byte[] bases=sequences.get(g.seqid);
		require(bases!=null, "GFF seqid absent from FASTA: "+g.seqid);
		require(g.start>=1 && g.stop>=g.start && g.stop<=bases.length, "Invalid 1-based inclusive GFF interval: "+coordinate(g));
	}
	/** Reference-only virtual intervals are recognized, never translated or remapped into primary stops. */
	private static boolean validateInterval(GffLine g, HashMap<String,byte[]> sequences, HashSet<String> circular){
		final byte[] bases=sequences.get(g.seqid);
		require(bases!=null, "GFF seqid absent from FASTA: "+g.seqid);
		if(g.stop<=bases.length){validateBounds(g, sequences); return false;}
		require(circular.contains(g.seqid) && g.start>=1 && g.start<=bases.length
			&& (long)g.stop-g.start+1<=bases.length,
			"Invalid circular-origin interval or missing validated region: "+coordinate(g));
		return true;
	}
	private static void addMask(Dataset data, GffLine g, int strand, HashMap<String,byte[]> sequences){
		final int n=sequences.get(g.seqid).length;
		assert(g.start>=1 && g.start<=n && g.stop>=g.start) : "Mask inputs must pass interval validation first";
		data.masks.add(new Mask(g.seqid, strand, g.start, Math.min(g.stop, n)));
		if(g.stop>n){data.masks.add(new Mask(g.seqid, strand, 1, g.stop-n));}
	}
	private static void writeContigs(String path, HashMap<String,byte[]> sequences, HashSet<String> circular){
		final ByteStreamWriter writer=new ByteStreamWriter(path, false, false, false); writer.start();
		try{
			writer.println("contig\tfasta_length\treference_circular");
			for(String id : new TreeSet<String>(sequences.keySet())){
				writer.println(id+"\t"+sequences.get(id).length+"\t"+circular.contains(id));
			}
		}finally{require(!writer.poisonAndWait(), "Failed writing FASTA contig inventory");}
	}
	private static boolean pseudo(GffLine g){
		final String value=attr(g, "pseudo"), type=attr(g, "gene_biotype");
		return (value!=null && !value.equals("false")) || attr(g, "pseudogene")!=null || "pseudogene".equals(type);
	}
	private static String attr(GffLine g, String key){return CdsTranslator.attribute(g.attributes, key);}
	private static String coordinate(GffLine g){return g.seqid+"\t"+g.strand+"\t"+g.start+"\t"+g.stop;}
	private static String stop(GffLine g){return g.seqid+"\t"+g.strand+"\t"+(g.strand==GffLine.MINUS ? g.start : g.stop);}
	private static void require(boolean condition, String message){if(!condition){throw new IllegalArgumentException(message);}}

	private static void writeRows(String path, Dataset data){
		final ByteStreamWriter writer=new ByteStreamWriter(path, false, false, false); writer.start();
		try{
			writer.println("schema\tinput_row\ttype\tcontig\tstrand\tstart\tend\tphase\ttable\tpartial5\tpartial3\tstratum\tattributes");
			for(Row row : data.rows){writeRow(writer, row);}
			for(Row row : data.maskRows){writeRow(writer, row);}
		}finally{require(!writer.poisonAndWait(), "Failed writing normalized GFF rows: "+path);}
	}
	private static void writeRow(ByteStreamWriter writer, Row row){
		final GffLine g=row.g;
		final ByteBuilder b=new ByteBuilder();
		b.append(SCHEMA).tab().append(row.ordinal).tab().append(g.type).tab().append(g.seqid).tab().append(g.strand).tab()
			.append(g.start).tab().append(g.stop).tab().append(g.phase).tab().append(row.table).tab().append(row.five ? 1 : 0).tab()
			.append(row.three ? 1 : 0).tab().append(row.reason).tab().append(g.attributes==null ? "." : g.attributes);
		writer.println(b);
	}

	private static void writeStops(String path, Dataset ref, Dataset before, Dataset after){
		final TreeSet<String> keys=new TreeSet<String>(); keys.addAll(ref.all.keySet()); keys.addAll(before.all.keySet()); keys.addAll(after.all.keySet());
		final ByteStreamWriter writer=new ByteStreamWriter(path, false, false, false); writer.start();
		try{
			writer.println("schema\tcontig\tstrand\tterminal\treference_raw_rows\treference_primary_rows\tbaseline_raw_rows\tbaseline_primary_rows\tcandidate_raw_rows\tcandidate_primary_rows\tchange\treference_rows\tbaseline_rows\tcandidate_rows");
			for(String key : keys){
				final ArrayList<Row> old=before.primary.get(key), next=after.primary.get(key);
				final String change=old==null ? (next==null ? "excluded_or_uncalled" : "gained") : next==null ? "lost"
					: coordinates(old).equals(coordinates(next)) ? "unchanged" : "same_stop_changed_starts";
				final ByteBuilder b=new ByteBuilder();
				b.append(SCHEMA).tab().append(key).tab().append(size(ref.all.get(key))).tab().append(size(ref.primary.get(key))).tab()
					.append(size(before.all.get(key))).tab().append(size(old)).tab().append(size(after.all.get(key))).tab().append(size(next)).tab()
					.append(change).tab().append(rowIds(ref.all.get(key))).tab().append(rowIds(before.all.get(key))).tab().append(rowIds(after.all.get(key)));
				writer.println(b);
			}
		}finally{require(!writer.poisonAndWait(), "Failed writing paired stop inventory");}
	}
	private static HashSet<String> coordinates(ArrayList<Row> rows){
		final HashSet<String> set=new HashSet<String>(); for(Row row : rows){set.add(coordinate(row.g));} return set;
	}
	private static int size(ArrayList<Row> rows){return rows==null ? 0 : rows.size();}
	private static String rowIds(ArrayList<Row> rows){
		if(rows==null){return ".";}
		final ByteBuilder b=new ByteBuilder();
		for(Row row : rows){if(b.length()>0){b.append(',');} b.append(row.ordinal);}
		return b.toString();
	}
	private static void writeSummary(ByteStreamWriter writer, String arm, Dataset data, Dataset ref){
		int primaryRows=0, five=0, matched=0;
		final TreeMap<String,Count> strata=new TreeMap<String,Count>();
		final TreeMap<String,Count> tables=new TreeMap<String,Count>();
		final HashSet<String> coords=new HashSet<String>();
		for(Row row : data.rows){
			coords.add(coordinate(row.g));
			Count count=strata.get(row.reason); if(count==null){count=new Count(); strata.put(row.reason, count);} count.value++;
			final String table=(attr(row.g, "transl_table")==null ? "fallback_" : "declared_")+row.table;
			Count tc=tables.get(table); if(tc==null){tc=new Count(); tables.put(table, tc);} tc.value++;
			if(row.primary()){primaryRows++; if(row.five){five++;}}
		}
		for(String key : data.primary.keySet()){if(ref.primary.containsKey(key)){matched++;}}
		assert(matched<=data.primary.size() && matched<=ref.primary.size() && primaryRows<=data.rawCds)
			: "Distinct-stop intersection cannot exceed either denominator and raw rows partition into primary/excluded";
		metric(writer, arm, "raw_cds_rows", data.rawCds); metric(writer, arm, "raw_coordinate_keys", coords.size());
		metric(writer, arm, "duplicate_coordinate_rows", data.rawCds-coords.size());
		metric(writer, arm, "raw_terminal_keys", data.all.size());
		metric(writer, arm, "excluded_rows", data.rawCds-primaryRows); metric(writer, arm, "primary_rows", primaryRows);
		metric(writer, arm, "primary_5prime_partial_rows", five); metric(writer, arm, "primary_stops", data.primary.size());
		metric(writer, arm, "reference_primary_stops", ref.primary.size()); metric(writer, arm, "matched_stops", matched);
		metric(writer, arm, "unmatched_primary_stops", data.primary.size()-matched); metric(writer, arm, "missed_reference_stops", ref.primary.size()-matched);
		int gross=0; for(String key : data.all.keySet()){if(!ref.primary.containsKey(key)){gross++;}}
		metric(writer, arm, "raw_terminal_keys_outside_primary_reference", gross); metric(writer, arm, "reference_mask_intervals", ref.masks.size());
		metric(writer, arm, "reference_pseudo_mask_rows", ref.maskRows.size());
		metric(writer, arm, "stop_reference_precision", ratio(matched, data.primary.size()));
		metric(writer, arm, "stop_reference_recall", ratio(matched, ref.primary.size()));
		for(String stratum : strata.keySet()){metric(writer, arm, "rows."+stratum, strata.get(stratum).value);}
		for(String table : tables.keySet()){metric(writer, arm, "table_rows."+table, tables.get(table).value);}
	}
	private static void writeDeltaSummary(ByteStreamWriter writer, Dataset ref, Dataset before, Dataset after){
		int gained=0, lost=0, gainedRef=0, lostRef=0, shared=0, starts=0;
		for(String key : before.primary.keySet()){
			if(!after.primary.containsKey(key)){lost++; if(ref.primary.containsKey(key)){lostRef++;}}
			else{shared++; if(!coordinates(before.primary.get(key)).equals(coordinates(after.primary.get(key)))){starts++;}}
		}
		for(String key : after.primary.keySet()){
			if(!before.primary.containsKey(key)){gained++; if(ref.primary.containsKey(key)){gainedRef++;}}
		}
		assert(shared+lost==before.primary.size() && shared+gained==after.primary.size())
			: "Paired stop partitions must conserve both arms, including start-only changes";
		metric(writer, "delta", "shared_primary_stops", shared); metric(writer, "delta", "gained_primary_stops", gained);
		metric(writer, "delta", "lost_primary_stops", lost); metric(writer, "delta", "gained_reference_stops", gainedRef);
		metric(writer, "delta", "lost_reference_stops", lostRef); metric(writer, "delta", "same_stop_changed_starts", starts);
	}
	private static String ratio(int n, int d){return d==0 ? "NA" : Double.toString((double)n/d);}
	private static void metric(ByteStreamWriter writer, String arm, String name, int value){metric(writer, arm, name, Integer.toString(value));}
	private static void metric(ByteStreamWriter writer, String arm, String name, String value){writer.println(SCHEMA+"\t"+arm+"\t"+name+"\t"+value);}

	public static final class Dataset {
		public Dataset(boolean reference_){reference=reference_;}
		public void index(){
			all.clear(); primary.clear();
			for(Row row : rows){add(all, row); if(row.primary()){add(primary, row);}}
			assert(primary.size()<=all.size()) : "Primary stop sites must be drawn from the raw input endpoint inventory";
		}
		private static void add(TreeMap<String,ArrayList<Row>> map, Row row){
			ArrayList<Row> list=map.get(row.stop); if(list==null){list=new ArrayList<Row>(); map.put(row.stop, list);} list.add(row);
		}
		final boolean reference;
		public int rawCds;
		public final ArrayList<Row> rows=new ArrayList<Row>(), maskRows=new ArrayList<Row>();
		public final ArrayList<Mask> masks=new ArrayList<Mask>();
		final HashSet<String> circular=new HashSet<String>();
		public final TreeMap<String,ArrayList<Row>> all=new TreeMap<String,ArrayList<Row>>(), primary=new TreeMap<String,ArrayList<Row>>();
	}
	public static final class Row {
		public Row(GffLine g_, int ordinal_){g=g_; ordinal=ordinal_; stop=stop(g);}
		public boolean primary(){return "primary".equals(reason);}
		public final GffLine g;
		final int ordinal; final String stop;
		GffLine parent, circularParent; Group group;
		public String reason;
		public int table=-1;
		public boolean five, three;
		boolean circular;
	}
	private static final class Group {final ArrayList<Row> rows=new ArrayList<Row>();}
	private static final class Count {int value;}
	private static final class Mask {
		Mask(String contig_, int strand_, int low_, int high_){contig=contig_; strand=strand_; low=low_; high=high_;}
		boolean contains(GffLine g){
			final int end=g.strand==GffLine.MINUS ? g.start : g.stop;
			return contig.equals(g.seqid) && (strand>1 || strand<0 || strand==g.strand) && end>=low && end<=high;
		}
		final String contig; final int strand, low, high;
	}
	static final String SCHEMA="stops-v1";
	private ProteinStopEvaluation(){}
}
