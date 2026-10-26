package test.assemble;

import java.io.PrintStream;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Locale;
import java.util.Map;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Reduces a complete matched Hollow suffix panel without dropping adverse genomes.
 * Nonexact suffixes are fixed-coordinate disagreements, not inferred edit counts.
 * @author Brian Bushnell, Fischl */
public final class HollowPanelSummary {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Uses the standard BBTools preparser, including config files. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, HollowPanelSummary.class, false);
		try{new HollowPanelSummary(pp.args, pp.outstream).process();}
		finally{Shared.closeStream(pp.outstream);}
	}

	/** Requires fresh output files and an explicit audit root. */
	private HollowPanelSummary(final String[] args, final PrintStream log_){
		log=log_;
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(String arg : args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(a.equals("root")){root=b;}
			else if(a.equals("details")){pooled=b;}
			else if(a.equals("cases")){cases=b;}
			else if(a.equals("genomes")){expectedGenomes=Integer.parseInt(b);}
			else if(a.equals("test")){test=Parse.parseBoolean(b);}
			else if(a.equals("pairmanifest")){pairManifest=Parse.parseBoolean(b);}
			else if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		input=parser.in1;
		output=parser.out1;
		require(!parser.overwrite && !parser.append, "Panel summaries must be fresh, never appended.");
		require(test || (input!=null && (root!=null || pairManifest) && output!=null && pooled!=null && expectedGenomes>0),
				"Specify in=manifest, root=audit_root, out=paired.tsv, details=pooled.tsv and genomes.");
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Loads all required arms from the fixed manifest; missing arms are fatal. */
	private void process(){
		if(test){selfTest(); return;}
		if(pairManifest){processPairs(); return;}
		final long[][][][] sums=new long[SETTINGS.length][2][2][8];
		final long[][][] changes=new long[SETTINGS.length][2][3];
		final HashSet<String> genomes=new HashSet<String>();
		final LineParser1 lp=new LineParser1('\t');
		final ByteBuilder row=new ByteBuilder();
		final ByteFile manifest=ByteFile.makeByteFile(input, true);
		try(Report pairs=new Report(output); Report census=cases==null ? null : new Report(cases)){
			pairs.writer.print(header("genome\tsetting\tcap", false));
			if(census!=null){census.writer.print(caseHeader());}
			for(byte[] line=manifest.nextLine(); line!=null; line=manifest.nextLine()){
				if(line.length==0 || line[0]=='#'){continue;}
				lp.set(line);
				final String genome=lp.parseString(0);
				require(genome.matches("g[0-9]{4}") && genomes.add(genome), "Invalid or duplicate genome: "+genome);
				for(int s=0; s<SETTINGS.length; s++){
					for(int c=0; c<CAPS.length; c++){
						final long[][] values=new long[2][];
						for(int g=0; g<2; g++){
							final String arm=genome+"_"+SETTINGS[s]+"_"+GEOMETRIES[g]+CAPS[c];
							values[g]=readSummary(root+"/runs/"+arm+"/summary.tsv", arm, CAPS[c]);
							for(int i=0; i<8; i++){sums[s][c][g][i]+=values[g][i];}
						}
						require(values[0][0]==values[1][0], "Matched geometries must contain the same seeds: "+genome);
						final long delta=nonexact(values[1])-nonexact(values[0]);
						changes[s][c][delta>0 ? 0 : delta<0 ? 1 : 2]++;
						row.clear().append(genome).tab().append(SETTINGS[s]).tab().append(CAPS[c]);
						appendPair(row, values);
						pairs.writer.print(row.nl());
						if(census!=null){writeCases(genome, SETTINGS[s], CAPS[c], values, census.writer);}
					}
				}
			}
			require(genomes.size()==expectedGenomes, "Manifest genome count differs from declared complete panel: "+genomes.size());
		}finally{require(!manifest.close(), "Manifest read failed.");}
		try(Report totals=new Report(pooled)){
			totals.writer.print(header("setting\tcap", true));
			for(int s=0; s<SETTINGS.length; s++){
				for(int c=0; c<CAPS.length; c++){
					row.clear().append(SETTINGS[s]).tab().append(CAPS[c]);
					appendPair(row, sums[s][c]);
					for(long change : changes[s][c]){row.tab().append(change);}
					totals.writer.print(row.nl());
				}
			}
		}
		log.println("HOLLOW_PANEL_SUMMARY_PASS genomes="+genomes.size()+" arms="+(genomes.size()*28));
		if(cases!=null){
			for(int i=0; i<caseTotals.length; i++){
				log.println("HOLLOW_PANEL_CASE_COUNT "+CASE_LABELS[i]+"="+caseTotals[i]);
			}
		}
	}

	/** Compares explicit baseline/treatment directories without renaming or copying frozen audits. */
	private void processPairs(){
		final HashSet<String> genomes=new HashSet<String>();
		final long[][] sums=new long[2][8];
		final long[] changes=new long[3];
		final String[] names={"baseline", "treatment"};
		final String[] labels={"both_nonexact", "baseline_only_nonexact", "treatment_only_nonexact", "circular_compatible_only"};
		final LineParser1 lp=new LineParser1('\t');
		final ByteBuilder row=new ByteBuilder();
		int panelCap=-1;
		final ByteFile manifest=ByteFile.makeByteFile(input, true);
		try(Report pairs=new Report(output); Report census=cases==null ? null : new Report(cases)){
			pairs.writer.print(header("genome\tcap", false, names));
			if(census!=null){census.writer.print(caseHeader("genome\tcap", names));}
			final byte[] schema=manifest.nextLine();
			require(schema!=null && new String(schema, java.nio.charset.StandardCharsets.US_ASCII).equals(PAIR_SCHEMA),
					"Explicit pairs require the documented six-column manifest header.");
			for(byte[] line=manifest.nextLine(); line!=null; line=manifest.nextLine()){
				lp.set(line);
				require(lp.terms()==6, "Explicit pairs require six columns per genome.");
				final String genome=lp.parseString(0);
				final int cap=lp.parseInt(1);
				require(genome.matches("g[0-9]{4}") && genomes.add(genome), "Invalid or duplicate genome: "+genome);
				if(panelCap<0){panelCap=cap;}
				require(cap>0 && cap==panelCap, "One paired report must use one positive extension cap.");
				final String[] paths={lp.parseString(3), lp.parseString(5)};
				final long[][] values=new long[2][];
				for(int g=0; g<2; g++){
					values[g]=readSummary(paths[g]+"/summary.tsv", lp.parseString(2+g*2), cap);
					for(int i=0; i<8; i++){sums[g][i]+=values[g][i];}
				}
				require(values[0][0]==values[1][0], "Paired treatments must contain the same seeds: "+genome);
				final long delta=nonexact(values[1])-nonexact(values[0]);
				changes[delta>0 ? 0 : delta<0 ? 1 : 2]++;
				row.clear().append(genome).tab().append(cap);
				appendPair(row, values);
				pairs.writer.print(row.nl());
				if(census!=null){writeCasePair(genome+"\t"+cap, cap, values, paths, labels, census.writer);}
			}
			require(genomes.size()==expectedGenomes, "Explicit pair count differs from declared panel: "+genomes.size());
		}finally{require(!manifest.close(), "Pair manifest read failed.");}
		try(Report totals=new Report(pooled)){
			totals.writer.print(header("cap", true, names));
			row.clear().append(panelCap);
			appendPair(row, sums);
			for(long value : changes){row.tab().append(value);}
			totals.writer.print(row.nl());
		}
		log.println("HOLLOW_EXPLICIT_PAIR_PASS genomes="+genomes.size()+" arms="+(genomes.size()*2));
		if(cases!=null){
			for(int i=0; i<labels.length; i++){log.println("HOLLOW_PANEL_CASE_COUNT "+labels[i]+"="+caseTotals[i]);}
		}
	}

	/** Joins both adverse reports, retaining Hollow-only failures even when totals improve. */
	private void writeCases(final String genome, final String setting, final int cap,
			final long[][] summaries, final ByteStreamWriter writer){
		final String[] paths=new String[2];
		for(int g=0; g<2; g++){
			final String arm=genome+"_"+setting+"_"+GEOMETRIES[g]+cap;
			paths[g]=root+"/runs/"+arm;
		}
		writeCasePair(genome+"\t"+setting+"\t"+cap, cap, summaries, paths, CASE_LABELS, writer);
	}

	/** Shares detail validation and ID pairing between the fixed panel and explicit treatments. */
	private void writeCasePair(final String prefix, final int cap, final long[][] summaries,
			final String[] paths, final String[] labels, final ByteStreamWriter writer){
		assert(summaries.length==2 && paths.length==2 && labels.length==4) : "Case comparisons require two arms and four disjoint classifications.";
		final LinkedHashMap<String, AuditRow[]> rows=new LinkedHashMap<String, AuditRow[]>();
		for(int g=0; g<2; g++){readCases(paths[g]+"/suffixes.tsv", g, cap, summaries[g], rows);}
		final ByteBuilder row=new ByteBuilder();
		for(Map.Entry<String, AuditRow[]> entry : rows.entrySet()){
			final AuditRow[] pair=entry.getValue();
			caseTotals[caseIndex(pair[0], pair[1])]++;
			row.clear().append(prefix).tab().append(entry.getKey()).tab().append(labels[caseIndex(pair[0], pair[1])]);
			appendCase(row, pair[0]);
			appendCase(row, pair[1]);
			writer.print(row.nl());
		}
	}

	/** Validates each detail file against its accepted summary before returning its retained rows. */
	private static void readCases(final String path, final int geometry, final int cap,
			final long[] summary, final LinkedHashMap<String, AuditRow[]> rows){
		assert(geometry>=0 && geometry<2 && summary.length==8) : "Case rows must belong to one validated panel arm.";
		final ByteFile file=ByteFile.makeByteFile(path, true);
		final LineParser1 lp=new LineParser1('\t');
		long internal=0, boundary=0, circular=0, disagreements=0;
		try{
			final byte[] header=file.nextLine();
			require(header!=null && new String(header, java.nio.charset.StandardCharsets.US_ASCII).equals(CASE_SCHEMA),
					"Unexpected case schema: "+path);
			for(byte[] line=file.nextLine(); line!=null; line=file.nextLine()){
				lp.set(line);
				require(lp.terms()==9, "Malformed suffix detail: "+path);
				final String id=lp.parseString(0);
				final AuditRow value=new AuditRow(lp);
				value.validate(cap);
				AuditRow[] pair=rows.get(id);
				if(pair==null){pair=new AuditRow[2]; rows.put(id, pair);}
				require(pair[geometry]==null, "Duplicate suffix ID in "+path+": "+id);
				pair[geometry]=value;
				if(value.beyond>0){boundary++;}
				else{internal++;}
				if(!value.nonexact()){circular++;}
				disagreements+=value.linear;
			}
			require(internal==summary[4] && boundary==summary[5] && circular==summary[6] && disagreements==summary[7],
					"Detail rows disagree with summary counters: "+path);
		}finally{require(!file.close(), "Suffix detail read failed: "+path);}
	}

	/** Absence from an adverse report means exact OR unextended; it never supplies an extension length. */
	private static String caseLabel(final AuditRow solid, final AuditRow hollow){
		return CASE_LABELS[caseIndex(solid, hollow)];
	}

	/** Keeps output labels and primitive census counters on the same partition. */
	private static int caseIndex(final AuditRow solid, final AuditRow hollow){
		assert(solid!=null || hollow!=null) : "The census emits only IDs present in at least one adverse report.";
		final boolean a=solid!=null && solid.nonexact(), b=hollow!=null && hollow.nonexact();
		return a ? (b ? 0 : 1) : (b ? 2 : 3);
	}

	/** Writes the original audit fields; -1 explicitly marks unavailable, not zero, counts. */
	private static void appendCase(final ByteBuilder row, final AuditRow value){
		assert(row!=null) : "The census owns one reusable output buffer.";
		if(value==null){row.append("\t-1\t-1\t-1\t-1\t-1\tnot_reported\t.\t."); return;}
		row.tab().append(value.added).tab().append(value.linear).tab().append(value.beyond)
				.tab().append(value.wrapped).tab().append(value.first).tab().append(value.classification)
				.tab().append(value.observed).tab().append(value.truth);
	}

	/** Names the retained evidence explicitly so downstream tools cannot mistake it for MSDI. */
	private static ByteBuilder caseHeader(){
		return caseHeader("genome\tsetting\tcap", GEOMETRIES);
	}

	/** Uses neutral baseline/treatment names for explicit pairs, preserving legacy column names. */
	private static ByteBuilder caseHeader(final String prefix, final String[] names){
		assert(names.length==2) : "Each comparison emits two sets of suffix fields.";
		final ByteBuilder row=new ByteBuilder().append(prefix).append("\tid\tcomparison");
		for(String geometry : names){
			for(String field : CASE_FIELDS){row.tab().append(geometry).append('_').append(field);}
		}
		return row.nl();
	}

	/** Requires exact schema, label, one data row and count conservation. */
	private static long[] readSummary(final String path, final String arm, final int cap){
		final ByteFile file=ByteFile.makeByteFile(path, true);
		final LineParser1 lp=new LineParser1('\t');
		try{
			final byte[] first=file.nextLine(), line=file.nextLine();
			require(first!=null && line!=null && new String(first, java.nio.charset.StandardCharsets.US_ASCII).equals(SCHEMA),
					"Unexpected suffix-audit schema: "+path);
			lp.set(line);
			require(lp.terms()==9 && lp.termEquals(arm, 0), "Wrong arm or column count: "+path);
			final long[] counts=new long[8];
			for(int i=0; i<counts.length; i++){counts[i]=lp.parseLong(i+1);}
			validate(counts, cap);
			require(file.nextLine()==null, "Multiple data rows in a single-arm summary: "+path);
			return counts;
		}finally{require(!file.close(), "Summary read failed: "+path);}
	}

	/** Mirrors the audit partition: exact, internal disagreement, or boundary. */
	private static void validate(final long[] v, final int cap){
		require(v.length==8 && cap>0, "Expected eight audit counters and a positive extension cap.");
		for(long x : v){require(x>=0, "Audit counters cannot be negative.");}
		require(v[0]>=v[1] && v[2]>=v[1] && v[2]<=v[1]*cap, "Seed/extension/base counts disagree.");
		require(v[3]+v[4]+v[5]==v[1] && v[6]<=v[5] && v[7]>=v[4], "Audit classes do not conserve extended seeds.");
	}

	/** Counts unverified suffixes, including boundary disagreements, without asserting topology. */
	private static long nonexact(final long[] v){
		assert(v.length==8 && v[5]>=v[6]) : "Circular-compatible seeds are a subset of boundary seeds.";
		return v[4]+v[5]-v[6];
	}

	/** Emits both raw counter sets before derived deltas, preserving adverse evidence. */
	private static void appendPair(final ByteBuilder row, final long[][] pair){
		assert(pair.length==2) : "The benchmark has exactly solid and Hollow arms.";
		for(long[] v : pair){for(long x : v){row.tab().append(x);}}
		row.tab().append(nonexact(pair[0])).tab().append(nonexact(pair[1]))
				.tab().append(pair[1][2]-pair[0][2]).tab().append(nonexact(pair[1])-nonexact(pair[0]));
	}

	/** Builds column labels matching appendPair, with optional genome-level regression counts. */
	private static ByteBuilder header(final String prefix, final boolean pooled){
		return header(prefix, pooled, GEOMETRIES);
	}

	/** Builds paired summary labels for either legacy geometries or explicit treatments. */
	private static ByteBuilder header(final String prefix, final boolean pooled, final String[] names){
		assert(names.length==2) : "appendPair writes exactly two counter sets.";
		final ByteBuilder bb=new ByteBuilder().append(prefix);
		for(String geometry : names){for(String field : FIELDS){bb.tab().append(geometry).append('_').append(field);}}
		bb.tab().append(names[0]).append("_nonexact\t").append(names[1]).append("_nonexact\tadded_bases_delta\tnonexact_delta");
		if(pooled){bb.append("\tgenomes_more_nonexact\tgenomes_fewer_nonexact\tgenomes_equal_nonexact");}
		return bb.nl();
	}

	/** Checks the class partition and ensures malformed summaries fail loudly. */
	private static void selfTest(){
		final long[] v={10, 8, 80, 4, 2, 2, 1, 3};
		validate(v, 20);
		require(nonexact(v)==3, "Boundary compatibility must not erase internal disagreements.");
		v[3]++;
		boolean rejected=false;
		try{validate(v, 20);}catch(IllegalArgumentException expected){rejected=true;}
		require(rejected, "Invalid class totals were accepted.");
		System.out.println("HOLLOW_PANEL_SUMMARY_TEST_PASS checks=2");
		final LineParser1 lp=new LineParser1('\t');
		lp.set("id\t3\t1\t0\t0\t1\tinternal_disagreement\tACA\tAAA".getBytes(java.nio.charset.StandardCharsets.US_ASCII));
		final AuditRow bad=new AuditRow(lp);
		bad.validate(20);
		require(caseLabel(null, bad).equals("hollow_only_nonexact"), "A new Hollow failure must not disappear behind improved totals.");
		require(caseLabel(bad, bad).equals("both_nonexact"), "Shared failures must remain identifiable.");
		lp.set("id\t3\t0\t1\t0\t-1\tcircular_compatible\tAAA\tAAA".getBytes(java.nio.charset.StandardCharsets.US_ASCII));
		final AuditRow edge=new AuditRow(lp);
		edge.validate(20);
		require(caseLabel(edge, bad).equals("hollow_only_nonexact") && caseLabel(null, edge).equals("circular_compatible_only"),
				"Circular compatibility must remain distinct from a known internal disagreement.");
		final ByteBuilder missing=new ByteBuilder();
		appendCase(missing, null);
		require(missing.toString().contains("not_reported"), "Missing adverse rows must not invent a zero-length solid extension.");
		System.out.println("HOLLOW_PANEL_CASE_TEST_PASS checks=4");
	}

	/** Retains only a rare adverse row, never per-base objects or the full read corpus. */
	private static final class AuditRow {
		/** Copies fields that must survive the next native input line. */
		AuditRow(final LineParser1 lp){
			assert(lp.terms()==9) : "The suffix audit defines nine columns, including its origin-labelled ID.";
			added=lp.parseInt(1); linear=lp.parseInt(2); beyond=lp.parseInt(3);
			wrapped=lp.parseInt(4); first=lp.parseInt(5);
			classification=lp.parseString(6);
			observed=lp.parseByteArray(7); truth=lp.parseByteArray(8);
		}

		/** Mirrors the existing audit partition without silently accepting exact rows. */
		void validate(final int cap){
			require(added>0 && added<=cap && observed.length==added && truth.length==added,
					"Suffix strings and extension cap must agree.");
			require(beyond>=0 && beyond<=added && linear>=0 && linear<=added-beyond && wrapped>=0 && wrapped<=beyond,
					"Suffix disagreement counts exceed their linear or wrapped domains.");
			require(first>=-1 && first<added && ((linear+wrapped==0)==(first==-1)), "First disagreement must agree with mismatch counts.");
			final String expected=beyond==0 ? "internal_disagreement" :
					(linear+wrapped==0 ? "circular_compatible" : "boundary_disagreement");
			require(classification.equals(expected) && (beyond>0 || linear>0), "Unexpected adverse-row classification.");
		}

		/** Circular-compatible rows are retained as boundary evidence, not counted as nonexact. */
		boolean nonexact(){return linear+wrapped>0;}

		final int added, linear, beyond, wrapped, first;
		final String classification;
		final byte[] observed, truth;
	}

	/** Keeps validation active even when assertions are disabled. */
	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	/** Owns a non-daemon report writer, including exceptional exit. */
	private static final class Report implements AutoCloseable {
		Report(final String path){
			writer=new ByteStreamWriter(path, false, false, true);
			writer.start();
		}
		@Override
		public void close(){require(!writer.poisonAndWait(), "Panel report write failed.");}
		final ByteStreamWriter writer;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String input, output, pooled, root, cases;
	private int expectedGenomes=100;
	private boolean test=false;
	private boolean pairManifest=false;
	private final PrintStream log;
	private final long[] caseTotals=new long[4];
	private static final String[] CASE_LABELS={"both_nonexact", "solid_only_nonexact", "hollow_only_nonexact", "circular_compatible_only"};
	private static final String[] SETTINGS={"control", "readins", "readdel", "readmixed", "seedsub", "seedins", "seeddel"};
	private static final int[] CAPS={20, 100};
	private static final String[] GEOMETRIES={"solid", "hollow"};
	private static final String[] FIELDS={"sequences", "extended", "added_bases", "linear_exact",
			"internal_disagreement", "boundary", "circular_compatible", "linear_disagreeing_bases"};
	private static final String SCHEMA="arm\tsequences\textended_sequences\tadded_bases\tlinear_exact\tinternal_disagreement"
			+"\tboundary_sequences\tcircular_compatible\tlinear_disagreeing_bases";
	private static final String CASE_SCHEMA="id\tadded\tlinear_disagreements\tbeyond_reference\twrapped_disagreements"
			+"\tfirst_disagreement_offset\tclassification\tobserved_suffix\tcoordinate_suffix";
	private static final String[] CASE_FIELDS={"added", "linear_disagreements", "beyond_reference", "wrapped_disagreements",
			"first_disagreement_offset", "classification", "observed_suffix", "coordinate_suffix"};
	private static final String PAIR_SCHEMA="genome\tcap\tbaseline_arm\tbaseline_dir\ttreatment_arm\ttreatment_dir";
}
