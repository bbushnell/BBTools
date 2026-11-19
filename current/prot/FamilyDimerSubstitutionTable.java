package prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Locale;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicInteger;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import parse.Parse;
import parse.PreParser;
import parse.LineParser1;
import shared.Shared;
import structures.ByteBuilder;

/** Observed adjacent-column dimer substitutions and fixed-odds weakest-link cuts.
 * Structural input verification is explicit; existing input byte identity is inherited.
 * @author Yaoyao
 */
public final class FamilyDimerSubstitutionTable{

	private FamilyDimerSubstitutionTable(){}
	public static final int DIMERS=400;
	public static final String[] BANDS={"<40", "40-60", ">=60"};
	public static final String[] STUDY_BANDS={"<40", "40-60", "60-80", ">=80"};

	public static void main(String[] args){
		try{
			final PreParser pp=new PreParser(args, FamilyDimerSubstitutionTable.class, false);
			try{run(pp.args);}finally{Shared.closeStream(pp.outstream);}
		}catch(Throwable failure){failure.printStackTrace(); System.exit(1);}
	}

	private static void run(String[] args) throws Exception{
		if(args.length>0 && args[0].equals("cuts")){cutArguments(Arrays.copyOfRange(args,1,args.length));return;}
		final HashMap<String,String> options=new HashMap<String,String>();
		for(String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected key=value option: "+arg);
			final String key=pair[0].toLowerCase(Locale.ROOT);
			require(OPTIONS.contains(key) && options.put(key, pair[1])==null, "Unknown or repeated option: "+key);
		}
		require(!options.containsKey("schema") || options.get("schema").equals("family_preparation_v1"),
			"Optional input schema must be family_preparation_v1");
		require(required(options, "verification").equalsIgnoreCase("structural"),
			"This reviewed recipe requires explicit verification=structural; no fingerprints are recomputed");
		final int threads=Math.toIntExact(number(options, "t"));
		require(threads>0 && threads<=64, "Use 1..64 workers for this bounded counting recipe");
		final int expectedFamilies=Math.toIntExact(number(options, "expectedfamilies"));
		final long expectedMembers=number(options, "expectedmembers");
		final Path out=Paths.get(required(options, "out")).toAbsolutePath().normalize();
		final IdentityBands customBands=options.containsKey("bands") ? IdentityBands.parse(options.get("bands")) : null;
		final Path control=customBands==null ? Paths.get(required(options, "monomercontrol")).toAbsolutePath().normalize() : null;
		require((control==null || Files.isRegularFile(control)) && Files.isDirectory(out.getParent()) && !Files.exists(out),
			"Use a readable accepted monomer control when required and a fresh output directory under an existing parent");
		Files.createDirectory(out);
		Shared.AMINO_IN=true; Shared.setThreads(threads);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=false;
		final long start=System.nanoTime();
		final FamilyCalibration.MemberInputs inputs=new FamilyCalibration.MemberInputs(
			Paths.get(required(options, "families")), required(options, "familiessha80"),
			Paths.get(required(options, "chunks")), required(options, "chunkssha80"),
			Paths.get(required(options, "consensus")), required(options, "consensussha80"),
			expectedFamilies, expectedMembers, Math.toIntExact(number(options, "buckets")),
			Math.toIntExact(number(options, "batches")), Math.toIntExact(number(options, "maxreferenceaa")),
			MemberLengthRange.from(options), FamilyCalibration.VerificationMode.STRUCTURAL);
		final long maxMembers=number(options, "maxbucketmembers"), maxBytes=number(options, "maxbucketbytes");
		final boolean study=options.containsKey("identitystudy") && Parse.parseBoolean(options.get("identitystudy"));
		require(!study || customBands==null, "identitystudy and custom bands are separate count modes");
		final Counts[] partial=new Counts[threads];
		for(int i=0; i<threads; i++){partial[i]=customBands==null ? new Counts(study) : new Counts(customBands);}
		final ExecutorService pool=Executors.newFixedThreadPool(threads);
		try{
			for(int bucket=0; bucket<inputs.buckets; bucket++){
				final FamilyCalibration.MemberInput[] families=inputs.loadBucket(bucket, maxMembers, maxBytes);
				final AtomicInteger next=new AtomicInteger();
				final ArrayList<Future<Void>> futures=new ArrayList<Future<Void>>(threads);
				for(int worker=0; worker<threads; worker++){futures.add(pool.submit(new Worker(families, next, partial[worker])));}
				for(Future<Void> future : futures){future.get();}
				if((bucket&15)==15 || bucket+1==inputs.buckets){System.err.println("Completed buckets: "+(bucket+1)+"/"+inputs.buckets);}
			}
		}finally{pool.shutdownNow();}
		inputs.verifyComplete();
		final Counts total=customBands==null ? new Counts(study) : new Counts(customBands);
		for(Counts counts : partial){total.add(counts);}
		total.verify();
		require(total.families==expectedFamilies && total.members==expectedMembers &&
			(customBands!=null || sum(total.alignments)==expectedMembers),
			"Every accepted member must be aligned exactly once, with every family accounted for");
		final Path monomer=out.resolve(customBands==null ? "monomer_control.tsv" : "monomer_counts.tsv");
		writeMonomer(monomer, customBands==null ? legacyCounts(total) : total);
		if(control!=null){
			require(Arrays.equals(Files.readAllBytes(monomer), Files.readAllBytes(control)),
				"Recomputed monomer control differs from accepted bytes; do not publish dimer maps");
		}
		if(study){writeIdentityStudy(out,total,args);return;}
		final long counted=System.nanoTime();
		final String[] labels=dimerLabels();
		final Table[] tables=new Table[total.bands.length];
		for(int b=0; b<total.bands.length; b++){tables[b]=new Table(total.ordered[b], DIMERS);}
		require(tables[0].total>0, "The raw identity<40 training band must contain observed dimers");
		final Hierarchy hierarchy=reduce(tables[0].odds, labels, 300, 200, 100);
		final long reduced=System.nanoTime();
		writeTables(out, total, tables, labels);
		writeHierarchy(out, hierarchy, labels, total.bands[0]);
		final StringBuilder summary=new StringBuilder("band\talignments\tmonomer_columns\tvalid_runs\tdimer_observations\tzero_endpoint_symbols\tzero_offdiagonal_cells\n");
		for(int b=0; b<total.bands.length; b++){
			int zeroSymbols=0, zeroPairs=0;
			for(long endpoint : tables[b].endpoints){if(endpoint==0){zeroSymbols++;}}
			for(int i=0; i<DIMERS; i++){for(int j=i+1; j<DIMERS; j++){if(tables[b].counts[i*DIMERS+j]==0){zeroPairs++;}}}
			summary.append(total.bands[b]).append('\t').append(total.alignments[b]).append('\t').append(total.columns[b])
				.append('\t').append(total.runs[b]).append('\t').append(total.dimers[b]).append('\t')
				.append(zeroSymbols).append('\t').append(zeroPairs).append('\n');
		}
		write(out.resolve("summary.tsv"), summary.toString());
		write(out.resolve("recipe.txt"), "schema=family_dimer_substitution_v1\nverification=structural\n"+
			"input_byte_identity=inherited_from_accepted_dataset_not_rehashed\n"+
			(customBands==null ? "control=byte_identical_original_monomer_table_not_proof_of_input_byte_identity\n" :
				"control=none_custom_identity_bands\n")+
			"alignment=StopAnchoredAminoLinear_QUERY_SUFFIX_traceback\n"+
			"adjacency=overlapping_consecutive_standard_m_columns_gap_or_unknown_breaks\n"+
			"training=raw_identity_lt40\nrule=maximize_min_original_odds\n"+
			"ties=exact_score_then_lexicographically_smallest_whole_postmerge_partition\n"+
			"serialization=numeric_member_order_comma_tokens_lexical_groups_slash_joined\n"+
			"count_seconds="+((counted-start)*1e-9)+"\nmerge_seconds="+((reduced-counted)*1e-9)+"\n"+String.join("\n", args)+"\n");
		write(out.resolve("PASS"), "FAMILY_DIMER_SUBSTITUTION_TABLE_PASS\n");
		System.out.println("FAMILY_DIMER_SUBSTITUTION_TABLE_PASS families="+total.families+" members="+total.members+
			" monomer_pairs="+sum(total.columns)+" dimer_pairs="+sum(total.dimers)+
			" control="+(customBands==null ? "byte_identical" : "not_required_custom_bands")+
			" hashes_computed=0 cuts=300,200,100 seconds="+((System.nanoTime()-start)*1e-9));
	}

	private static final class Worker implements Callable<Void>{
		Worker(FamilyCalibration.MemberInput[] families_, AtomicInteger next_, Counts counts_){families=families_; next=next_; counts=counts_;}
		@Override public Void call(){
			try{
				final FamilyCalibration.Buffers buffers=new FamilyCalibration.Buffers();
				for(int index=next.getAndIncrement(); index<families.length; index=next.getAndIncrement()){
					require(!Thread.currentThread().isInterrupted(), "Interrupted counting worker cannot publish a partial table");
					final FamilyCalibration.MemberInput family=families[index];
					final ProteinSequence consensus=new ProteinSequence(family.id, family.consensus);
					counts.families++; counts.members=Math.addExact(counts.members, family.members.length);
					for(byte[] raw : family.members){
						final ProteinSequence query=buffers.encode(raw);
						final AAAlignment alignment=StopAnchoredAminoLinear.align(query.enc, consensus.enc,
							StopAnchoredAminoLinear.Mode.QUERY_SUFFIX, true);
						count(query.enc, consensus.enc, alignment, counts);
					}
					families[index]=null;
				}
				return null;
			}catch(RuntimeException | Error failure){next.set(families.length); throw failure;}
		}
		final FamilyCalibration.MemberInput[] families;
		final AtomicInteger next;
		final Counts counts;
	}

	private static void writeMonomer(Path path, Counts counts){
		final ByteStreamWriter out=writer(path); final ByteBuilder row=new ByteBuilder();
		try{
			out.println("band\tquery_aa\tconsensus_aa\tcount\tquery_monomer\tconsensus_monomer\talignments\taligned_pairs");
			for(int b=0; b<counts.bands.length; b++){
				final long[] query=new long[20], consensus=new long[20];
				for(int q=0; q<20; q++){for(int t=0; t<20; t++){query[q]+=counts.monomer[b][20*q+t]; consensus[t]+=counts.monomer[b][20*q+t];}}
				for(int q=0; q<20; q++){
					for(int t=0; t<20; t++){
						row.clear().append(counts.bands[b]).tab().append((char)AminoAcid.numberToAcid[q]).tab()
							.append((char)AminoAcid.numberToAcid[t]).tab().append(counts.monomer[b][20*q+t]).tab()
							.append(query[q]).tab().append(consensus[t]).tab().append(counts.alignments[b]).tab().append(counts.columns[b]);
						out.println(row);
					}
				}
			}
		}finally{require(!out.poisonAndWait(), "Monomer control output failed");}
	}

	private static void writeTables(Path root, Counts counts, Table[] tables, String[] labels){
		final ByteStreamWriter dictionary=writer(root.resolve("dictionary.tsv"));
		try{dictionary.println("dimer_id\tdimer"); for(int i=0; i<labels.length; i++){dictionary.println(i+"\t"+labels[i]);}}
		finally{require(!dictionary.poisonAndWait(), "Dictionary output failed");}
		final ByteStreamWriter ordered=writer(root.resolve("ordered_counts.tsv"));
		final ByteBuilder row=new ByteBuilder();
		try{
			ordered.println("band\tquery_dimer\tconsensus_dimer\tcount");
			for(int b=0; b<counts.bands.length; b++){
				for(int i=0; i<DIMERS; i++){
					for(int j=0; j<DIMERS; j++){
						row.clear().append(counts.bands[b]).tab().append(labels[i]).tab()
							.append(labels[j]).tab().append(counts.ordered[b][i*DIMERS+j]);
						ordered.println(row);
					}
				}
			}
		}finally{require(!ordered.poisonAndWait(), "Ordered dimer output failed");}
		final ByteStreamWriter symmetric=writer(root.resolve("symmetric_tables.tsv"));
		try{
			symmetric.println("band\tdimer_i\tdimer_j\tcount\tendpoint_i\tendpoint_j\tobserved_frequency\texpected_frequency\todds");
			for(int b=0; b<counts.bands.length; b++){
				final Table table=tables[b];
				for(int i=0; i<DIMERS; i++){
					for(int j=i; j<DIMERS; j++){
						final long value=table.counts[i*DIMERS+j];
						final double pi=table.total==0 ? 0 : table.endpoints[i]/(2.0*table.total);
						final double pj=table.total==0 ? 0 : table.endpoints[j]/(2.0*table.total);
						row.clear().append(counts.bands[b]).tab().append(labels[i]).tab().append(labels[j]).tab().append(value).tab()
							.append(table.endpoints[i]).tab().append(table.endpoints[j]).tab()
							.append(Double.toString(table.total==0 ? 0 : value/(double)table.total)).tab()
							.append(Double.toString(i==j ? pi*pi : 2*pi*pj)).tab().append(Double.toString(table.odds[i*DIMERS+j]));
						symmetric.println(row);
					}
				}
			}
		}finally{require(!symmetric.poisonAndWait(), "Symmetric dimer output failed");}
	}

	private static void writeHierarchy(Path root, Hierarchy hierarchy, String[] labels){
		writeHierarchy(root, hierarchy, labels, "<40");
	}
	private static void writeHierarchy(Path root, Hierarchy hierarchy, String[] labels, String band){
		write(root.resolve("merge_ladder.tsv"), ladderText(hierarchy,0));
		write(root.resolve("merge_groups.tsv"), groupsText(hierarchy,false, band));
		write(root.resolve("maps.tsv"), mapsText(hierarchy,labels,false));
	}

	private static String ladderText(Hierarchy hierarchy,int minimum){
		final StringBuilder ladder=new StringBuilder("symbols_after_merge\tleft\tright\tjoined\tweakest_odds\n");
		for(Merge merge : hierarchy.merges){if(merge.symbols>=minimum){ladder.append(merge.symbols).append('\t').append(merge.a).append('\t').append(merge.b)
			.append('\t').append(merge.joined).append('\t').append(merge.odds).append('\n');}
		}
		return ladder.toString();
	}
	private static String groupsText(Hierarchy hierarchy,boolean controlsOnly){
		return groupsText(hierarchy, controlsOnly, "<40");
	}
	private static String groupsText(Hierarchy hierarchy,boolean controlsOnly,String band){
		final StringBuilder groups=new StringBuilder("weighting\tband\trule\tsymbols\tgroups\n");
		for(Cut cut : hierarchy.cuts){
			if(controlsOnly && cut.symbols!=300 && cut.symbols!=200 && cut.symbols!=100){continue;}
			groups.append("raw\t").append(band).append("\tweakest_link_odds\t").append(cut.symbols)
				.append('\t').append(String.join("/", cut.groups)).append('\n');
		}
		return groups.toString();
	}
	private static String mapsText(Hierarchy hierarchy,String[] labels,boolean controlsOnly){
		final StringBuilder maps=new StringBuilder("symbols\tdimer_id\tdimer\tclass_id\n");
		for(Cut cut : hierarchy.cuts){
			if(controlsOnly && cut.symbols!=300 && cut.symbols!=200 && cut.symbols!=100){continue;}
			for(int i=0; i<cut.map.length; i++){maps.append(cut.symbols).append('\t').append(i).append('\t').append(labels[i]).append('\t').append(cut.map[i]).append('\n');}
		}
		return maps.toString();
	}

	/** Reuse accepted original odds; the existing300 merges and three maps must reproduce exactly. */
	private static void cutArguments(String[] args)throws Exception{
		final HashMap<String,String> options=new HashMap<String,String>();
		for(String arg:args){
			final String[] pair=Parse.splitOnFirst(arg,'=');
			require(pair.length==2,"Cut export requires key=value arguments");
			final String key=pair[0].toLowerCase(Locale.ROOT);
			require(Arrays.asList("in","ladder","groups","maps","out","cuts","t").contains(key) && options.put(key,pair[1])==null,
				"Unknown or duplicate cut-export argument: "+key);
		}
		require(!options.containsKey("t") || number(options,"t")==1,"Cut export uses one thread and no member alignment");
		final String[] values=required(options,"cuts").split(",",-1);final int[] targets=new int[values.length];
		for(int i=0;i<targets.length;i++){targets[i]=Integer.parseInt(values[i]);}
		exportCuts(Paths.get(required(options,"in")),Paths.get(required(options,"ladder")),Paths.get(required(options,"groups")),
			Paths.get(required(options,"maps")),Paths.get(required(options,"out")),targets);
	}

	/** Only small controls are read twice; the full odds table is streamed once without hashing. */
	public static void exportCuts(Path input,Path ladder,Path groups,Path maps,Path out,int... targets)throws Exception{
		require(targets!=null && targets.length>0,"At least one requested hierarchy cut is required");
		final int[] all=new int[targets.length+3];all[0]=300;all[1]=200;all[2]=100;
		for(int i=0;i<targets.length;i++){
			require(targets[i]>0 && targets[i]<=DIMERS,"Requested cut must retain1..400 groups");
			for(int j=0;j<i;j++){require(targets[i]!=targets[j],"Duplicate requested cut");}
			all[i+3]=targets[i];
		}
		require(out.isAbsolute() && Files.isDirectory(out.getParent()) && !Files.exists(out) && !Files.isSymbolicLink(out),
			"Cut export requires a fresh absolute output directory");
		for(Path path:new Path[]{input,ladder,groups,maps}){require(path.isAbsolute() && Files.isRegularFile(path),"Require explicit readable cut input: "+path);}
		final java.nio.file.attribute.BasicFileAttributes before=Files.readAttributes(input,java.nio.file.attribute.BasicFileAttributes.class);
		final byte[] oldLadder=Files.readAllBytes(ladder),oldGroups=Files.readAllBytes(groups),oldMaps=Files.readAllBytes(maps);
		Shared.setThreads(1);ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=false;
		final String[] labels=dimerLabels();final double[] odds=readOdds(input,labels);
		final Hierarchy hierarchy=reduce(odds,labels,all);
		require(Arrays.equals(oldLadder,ladderText(hierarchy,100).getBytes(StandardCharsets.UTF_8)),"Accepted300-merge ladder differs; do not publish new cuts");
		require(Arrays.equals(oldGroups,groupsText(hierarchy,true).getBytes(StandardCharsets.UTF_8)),"Accepted300/200/100 groups differ");
		require(Arrays.equals(oldMaps,mapsText(hierarchy,labels,true).getBytes(StandardCharsets.UTF_8)),"Accepted300/200/100 maps differ");
		final java.nio.file.attribute.BasicFileAttributes after=Files.readAttributes(input,java.nio.file.attribute.BasicFileAttributes.class);
		require(before.size()==after.size() && before.lastModifiedTime().equals(after.lastModifiedTime()) &&
			java.util.Objects.equals(before.fileKey(),after.fileKey()),"Retained odds input changed during cut export");
		require(Arrays.equals(oldLadder,Files.readAllBytes(ladder)) && Arrays.equals(oldGroups,Files.readAllBytes(groups)) &&
			Arrays.equals(oldMaps,Files.readAllBytes(maps)),"Accepted hierarchy controls changed during export");
		Files.createDirectory(out);writeHierarchy(out,hierarchy,labels);
		write(out.resolve("recipe.txt"),"schema=family_dimer_cuts_v1\ntraining=raw_identity_lt40\nrule=maximize_min_original_odds\n"+
			"input_byte_identity=inherited_not_rehashed\ncontrol=exact300_merge_prefix_and300_200_100_maps\n"+
			"in="+input+"\nladder="+ladder+"\ngroups="+groups+"\nmaps="+maps+"\ncuts="+Arrays.toString(targets)+"\n");
		write(out.resolve("PASS"),"FAMILY_DIMER_CUTS_PASS\n");
		System.out.println("FAMILY_DIMER_CUTS_PASS merges="+hierarchy.merges.size()+" cuts="+hierarchy.cuts.size()+" hashes_computed=0");
	}

	private static double[] readOdds(Path input,String[] labels){
		final double[] odds=new double[DIMERS*DIMERS];
		final LineParser1 parser=new LineParser1('\t');
		final ByteFile reader=ByteFile.makeByteFile1(input.toString(),false);
		try{
			byte[] line=reader.nextLine();
			require(line!=null && new String(line,StandardCharsets.US_ASCII).equals(
				"band\tdimer_i\tdimer_j\tcount\tendpoint_i\tendpoint_j\tobserved_frequency\texpected_frequency\todds"),"Retained symmetric table schema differs");
			for(String band:BANDS){for(int i=0;i<DIMERS;i++){for(int j=i;j<DIMERS;j++){
				line=reader.nextLine();require(line!=null,"Truncated retained symmetric table");
				parser.set(line);
				require(parser.terms()==9 && parser.termEquals(band,0) && parser.termEquals(labels[i],1) && parser.termEquals(labels[j],2),
					"Retained odds rows must cover each band/upper-triangle coordinate once in dictionary order");
				// Exact Double.toString round-trip preserves merge scores and ties; the fast decimal parser need not.
				final double value=Double.parseDouble(parser.parseString(8));
				require(Double.isFinite(value) && value>=0,"Retained odds must be finite and nonnegative");
				if(band.equals("<40")){odds[i*DIMERS+j]=odds[j*DIMERS+i]=value;}
			}}}
			require(reader.nextLine()==null,"Unexpected extra retained odds rows");
		}finally{require(!reader.close(),"Retained odds stream failed");}
		return odds;
	}

	private static ByteStreamWriter writer(Path path){
		final ByteStreamWriter writer=new ByteStreamWriter(path.toString(), false, false, false);
		writer.start(); return writer;
	}
	private static void write(Path path, String text){
		try{Files.write(path, text.getBytes(StandardCharsets.UTF_8), java.nio.file.StandardOpenOption.CREATE_NEW);}
		catch(java.io.IOException failure){throw new RuntimeException("Output failed: "+path, failure);}
	}
	private static String required(HashMap<String,String> options, String key){
		final String value=options.get(key); require(value!=null && !value.isEmpty(), "Missing option: "+key); return value;
	}
	private static long number(HashMap<String,String> options, String key){
		final long value=Long.parseLong(required(options, key)); require(value>=0, "Negative count: "+key); return value;
	}
	private static final ArrayList<String> OPTIONS=new ArrayList<String>(Arrays.asList(
		"families", "familiessha80", "chunks", "chunkssha80", "consensus", "consensussha80", "expectedfamilies", "expectedmembers",
		"buckets", "batches", "maxreferenceaa", "maxbucketmembers", "maxbucketbytes", "minmemberaa", "maxmemberaa",
		"verification", "monomercontrol", "out", "t", "schema", "identitystudy", "bands"));

	/** Existing alignment identity includes gap columns; do not use dimer identity. */
	public static int band(double identity){
		require(Double.isFinite(identity) && identity>=0 && identity<=100, "Invalid alignment identity: "+identity);
		return identity<40 ? 0 : identity<60 ? 1 : 2;
	}


	/** Split the original high band without changing the identity definition. */
	public static int studyBand(double identity){
		final int old=band(identity);
		return old<2 ? old : identity<80 ? 2 : 3;
	}

	/** Collapse the extended bins for the byte-exact original monomer control. */
	public static Counts legacyCounts(Counts counts){
		require(counts!=null,"Counts are required for control projection");
		require(counts.customBands==null,"Overlapping custom bands cannot be projected to the three legacy bins");
		if(!counts.study){return counts;}
		final Counts out=new Counts();out.families=counts.families;out.members=counts.members;
		for(int b=0;b<4;b++){
			final int old=Math.min(b,2);
			out.columns[old]=Math.addExact(out.columns[old],counts.columns[b]);
			out.runs[old]=Math.addExact(out.runs[old],counts.runs[b]);
			out.dimers[old]=Math.addExact(out.dimers[old],counts.dimers[b]);
			out.alignments[old]=Math.addExact(out.alignments[old],counts.alignments[b]);
			for(int i=0;i<DIMERS*DIMERS;i++){out.ordered[old][i]=Math.addExact(out.ordered[old][i],counts.ordered[b][i]);}
			for(int i=0;i<400;i++){out.monomer[old][i]=Math.addExact(out.monomer[old][i],counts.monomer[b][i]);}
		}
		out.verify();return out;
	}

	/** Cumulative matrices include each observation exactly once, with strict upper cutoffs. */
	public static long[] cumulative(long[][] bands,int count){
		require(bands!=null && count>0 && count<=bands.length,"Cumulative cutoff must fit the disjoint bands");
		final long[] out=new long[bands[0].length];
		for(int b=0;b<count;b++){
			require(bands[b].length==out.length,"All identity bands use one complete dictionary");
			for(int i=0;i<out.length;i++){require(bands[b][i]>=0,"Observation counts cannot be negative");out[i]=Math.addExact(out[i],bands[b][i]);}
		}
		return out;
	}

	/** One alignment population, four cumulative cutoffs, three requested encodings. */
	public static void writeIdentityStudy(Path out,Counts counts,String[] args){
		require(counts!=null && counts.study,"Identity study requires the four disjoint identity bands");counts.verify();
		final String[] names={"lt40","lt60","lt80","all"};
		final String[] monomers=new String[20];
		for(int i=0;i<20;i++){monomers[i]=String.valueOf((char)AminoAcid.numberToAcid[i]);}
		final String[] dimers=dimerLabels();
		final StringBuilder groups=new StringBuilder("cutoff\tencoding\tsymbols\tgroups\n");
		final StringBuilder summary=new StringBuilder("cutoff\talignments\tmonomer_columns\tvalid_runs\tdimer_observations\n");
		long aligned=0,columns=0,runs=0,pairs=0;
		for(int b=0;b<4;b++){
			aligned+=counts.alignments[b];columns+=counts.columns[b];runs+=counts.runs[b];pairs+=counts.dimers[b];
			summary.append(names[b]).append('\t').append(aligned).append('\t').append(columns).append('\t').append(runs).append('\t').append(pairs).append('\n');
			final Table mono=new Table(cumulative(counts.monomer,b+1),20),di=new Table(cumulative(counts.ordered,b+1),400);
			require(mono.total>0 && di.total>0,"Every requested cutoff needs observed substitution evidence");
			final Hierarchy mh=reduce(mono.odds,monomers,"",16,12),dh=reduce(di.odds,dimers,120);
			for(Hierarchy h:new Hierarchy[]{mh,dh}){
				for(Cut cut:h.cuts){
					final String encoding=cut.symbols==120 ? "d120K6_nonoverlap" : cut.symbols==16 ? "16K5" : "12K6";
					groups.append(names[b]).append('\t').append(encoding).append('\t').append(cut.symbols).append('\t').append(String.join("/",cut.groups)).append('\n');
					write(out.resolve(names[b]+"_"+encoding+".groups.txt"),String.join("/",cut.groups)+"\n");
				}
			}
			write(out.resolve(names[b]+"_monomer_ladder.tsv"),ladderText(mh,0));
			write(out.resolve(names[b]+"_dimer_ladder.tsv"),ladderText(dh,0));
			write(out.resolve(names[b]+"_dimer_maps.tsv"),mapsText(dh,dimers,false));
		}
		write(out.resolve("merge_groups.tsv"),groups.toString());write(out.resolve("summary.tsv"),summary.toString());
		final ByteStreamWriter ordered=writer(out.resolve("study_ordered_counts.tsv.gz"));final ByteBuilder row=new ByteBuilder();
		try{
			ordered.println("band\tkind\tquery_symbol\tconsensus_symbol\tcount");
			for(int b=0;b<4;b++){
				for(int kind=0;kind<2;kind++){
					final String[] labels=kind==0 ? monomers : dimers;final long[] values=kind==0 ? counts.monomer[b] : counts.ordered[b];
					for(int i=0;i<labels.length;i++){for(int j=0;j<labels.length;j++){
						row.clear().append(STUDY_BANDS[b]).tab().append(kind==0 ? "monomer" : "dimer").tab().append(labels[i]).tab().append(labels[j]).tab().append(values[i*labels.length+j]);ordered.println(row);
					}}
				}
			}
		}finally{require(!ordered.poisonAndWait(),"Study ordered count output failed");}
		write(out.resolve("recipe.txt"),"schema=family_identity_cutoff_study_v1\ncutoffs=strict_lt40,strict_lt60,strict_lt80,all\n"+
			"population=cap150_member_to_consensus_not_all_pairs\nweighting=raw\nrule=maximize_min_original_odds\n"+
			"alignment=StopAnchoredAminoLinear_QUERY_SUFFIX_traceback\nidentity=AAAlignment.pident_including_gaps\n"+
			"dimer_observations=overlapping_adjacent_standard_columns\ndimer_encoding=nonoverlapping_K6_pairs_1_3_5\n"+String.join("\n",args)+"\n");
		write(out.resolve("PASS"),"FAMILY_IDENTITY_CUTOFF_STUDY_PASS\n");
		System.out.println("FAMILY_IDENTITY_CUTOFF_STUDY_PASS families="+counts.families+" members="+counts.members+" cutoffs=4 alphabets=12 monomer_control=byte_identical");
	}

	/** Every emitted dimer consumes two consecutive standard m columns in this trace.
	 * See AAAlignment.match and FamilySubstitutionPairTable.walk for coordinates.
	 */
	public static void count(byte[] query, byte[] target, AAAlignment alignment, Counts counts){
		require(query!=null && target!=null && alignment!=null && counts!=null && alignment.match!=null,
			"Dimer counting requires encoded sequences and a recorded alignment path");
		require(alignment.qStart>=0 && alignment.qStop<query.length && alignment.tStart>=0
			&& alignment.tStop<target.length && alignment.match.length==alignment.length,
			"Trace coordinates and length must refer to the supplied sequences");
		final int mask=counts.customBands==null ? 1<<(counts.study ? studyBand(alignment.pident()) : band(alignment.pident())) :
			counts.customBands.mask(alignment.pident());
		if(mask==0){return;}
		int q=alignment.qStart, t=alignment.tStart, previousQuery=-1, previousTarget=-1;
		boolean adjacent=false;
		for(byte op : alignment.match){
			if(op=='m'){
				assert(q<=alignment.qStop && t<=alignment.tStop) : "m consumes both sequence coordinates";
				final byte qa=query[q++], ta=target[t++];
				if(Blosum62.isStandard(qa) && Blosum62.isStandard(ta)){
					for(int band=0; band<counts.bands.length; band++){
						if((mask&(1<<band))!=0){
							counts.monomer[band][20*qa+ta]++;
							counts.columns[band]++;
							if(adjacent){
								final int queryDimer=20*previousQuery+qa, targetDimer=20*previousTarget+ta;
								counts.ordered[band][DIMERS*queryDimer+targetDimer]++;
								counts.dimers[band]++;
							}else{counts.runs[band]++;}
						}
					}
					previousQuery=qa; previousTarget=ta; adjacent=true;
				}else{adjacent=false;}
			}else if(op=='I'){q++; adjacent=false;}
			else if(op=='D'){t++; adjacent=false;}
			else{throw new IllegalArgumentException("Unknown alignment operation: "+(char)op);}
		}
		require(q==alignment.qStop+1 && t==alignment.tStop+1,
			"Trace must consume its declared endpoints; clipped prefixes cannot contribute dimers");
		for(int band=0; band<counts.bands.length; band++){
			if((mask&(1<<band))!=0){counts.alignments[band]++;}
		}
	}

	/** Reuse one instance per worker across buckets; never allocate one per family. */
	public static final class Counts{
		public Counts(){this(false, null);}
		public Counts(boolean study_){
			this(study_, null);
		}
		public Counts(IdentityBands customBands_){
			this(false, customBands_);
		}
		private Counts(boolean study_, IdentityBands customBands_){
			require(!study_ || customBands_==null, "Identity study cannot use overlapping custom bands");
			study=study_; customBands=customBands_;
			bands=customBands==null ? (study ? STUDY_BANDS : BANDS) : customBands.labels();
			ordered=new long[bands.length][DIMERS*DIMERS];monomer=new long[bands.length][400];
			columns=new long[bands.length];runs=new long[bands.length];
			dimers=new long[bands.length];alignments=new long[bands.length];
		}
		public final boolean study;
		public final IdentityBands customBands;
		public final String[] bands;
		public final long[][] ordered,monomer;
		public final long[] columns,runs,dimers,alignments;
		public long families, members;

		public void add(Counts other){
			require(other!=null && other!=this && other.study==study && other.customBands==customBands,
				"Worker reduction requires distinct counters with identical identity bins");
			families=Math.addExact(families, other.families);
			members=Math.addExact(members, other.members);
			for(int b=0; b<bands.length; b++){
				columns[b]=Math.addExact(columns[b], other.columns[b]);
				runs[b]=Math.addExact(runs[b], other.runs[b]);
				dimers[b]=Math.addExact(dimers[b], other.dimers[b]);
				alignments[b]=Math.addExact(alignments[b], other.alignments[b]);
				for(int i=0; i<ordered[b].length; i++){ordered[b][i]=Math.addExact(ordered[b][i], other.ordered[b][i]);}
				for(int i=0; i<monomer[b].length; i++){monomer[b][i]=Math.addExact(monomer[b][i], other.monomer[b][i]);}
			}
		}

		public void verify(){
			for(int b=0; b<bands.length; b++){
				require(columns[b]>=runs[b] && runs[b]>=0 && dimers[b]==columns[b]-runs[b],
					"Each positive eligible run of L columns emits L-1 dimers: band="+bands[b]);
				require(sum(ordered[b])==dimers[b] && sum(monomer[b])==columns[b],
					"Matrix totals must preserve the trace walker observations: band="+bands[b]);
			}
		}
	}

	/** Symmetrize endpoints, not the internal residue order of the dimer symbols. */
	public static final class Table{
		public Table(long[] ordered, int dimension){
			require(dimension>0 && dimension<=DIMERS && ordered!=null && ordered.length==dimension*dimension,
				"Ordered matrix dimensions must match the complete symbol dictionary");
			n=dimension; total=sum(ordered);
			counts=new long[n*n]; endpoints=new long[n]; odds=new double[n*n];
			for(int i=0; i<n; i++){
				for(int j=0; j<n; j++){
					final long value=ordered[i*n+j];
					endpoints[i]=Math.addExact(endpoints[i], value);
					endpoints[j]=Math.addExact(endpoints[j], value);
				}
			}
			require(sum(endpoints)==Math.multiplyExact(2L, total), "Every observation supplies two dimer endpoints");
			long symmetricTotal=0;
			for(int i=0; i<n; i++){
				final double pi=total==0 ? 0 : endpoints[i]/(2.0*total);
				for(int j=i; j<n; j++){
					final long count=i==j ? ordered[i*n+i] : Math.addExact(ordered[i*n+j], ordered[j*n+i]);
					counts[i*n+j]=count;
					symmetricTotal=Math.addExact(symmetricTotal, count);
					final double pj=total==0 ? 0 : endpoints[j]/(2.0*total);
					final double expected=(i==j ? pi*pi : 2*pi*pj);
					final double value=expected>0 ? (count/(double)total)/expected : 0;
					require(Double.isFinite(value) && value>=0, "Observed odds must be finite and nonnegative");
					odds[i*n+j]=odds[j*n+i]=value;
				}
			}
			require(symmetricTotal==total, "Symmetrization counts diagonals once and preserves all ordered observations");
		}
		public final int n;
		public final long total;
		/** Counts uses the upper triangle; odds is populated symmetrically. */
		public final long[] counts, endpoints;
		public final double[] odds;
	}

	public static String[] dimerLabels(){
		final String[] labels=new String[DIMERS];
		for(int a=0; a<20; a++){
			for(int b=0; b<20; b++){
				labels[20*a+b]=new String(new char[]{(char)AminoAcid.numberToAcid[a], (char)AminoAcid.numberToAcid[b]});
			}
		}
		return labels;
	}

	/** Complete-link updates preserve minima of the ORIGINAL odds, without pooling counts.
	 * Ties compare the entire sorted post-merge partition, including delimiters.
	 */
	public static Hierarchy reduce(double[] original, String[] labels, int... targets){
		return reduce(original, labels, ",", targets);
	}

	/** Empty member separator reproduces the original single-letter monomer tie key. */
	public static Hierarchy reduce(double[] original, String[] labels, String separator, int... targets){
		require(",".equals(separator) || "".equals(separator), "Only canonical dimer or legacy monomer serialization is supported");
		require(labels!=null && labels.length>0 && labels.length<=DIMERS && original!=null
			&& original.length==labels.length*labels.length && targets!=null && targets.length>0,
			"Reduction requires a complete square odds table, labels and requested cuts");
		final int n=labels.length;
		int minimum=n;
		for(int target : targets){require(target>=1 && target<=n, "Cut must fit the alphabet: "+target); minimum=Math.min(minimum, target);}
		final Group[] groups=new Group[n];
		final String[] mergedCache=new String[n*n];
		final double[] scores=original.clone();
		for(int i=0; i<n; i++){
			require(labels[i]!=null && labels[i].matches("[A-Z]+"), "Labels must be nonempty residue tokens without delimiters");
			require(!separator.isEmpty() || labels[i].length()==labels[0].length(),
				"Delimiter-free labels must have equal widths so distinct groups cannot serialize to the same text");
			for(int j=0; j<i; j++){require(!labels[i].equals(labels[j]), "Dictionary labels must be unique");}
			groups[i]=new Group(new int[]{i}, labels[i]);
			for(int j=0; j<n; j++){
				require(Double.isFinite(scores[i*n+j]) && scores[i*n+j]>=0 && scores[i*n+j]==scores[j*n+i],
					"Weakest-link odds must be finite, nonnegative and symmetric");
			}
		}
		final Hierarchy result=new Hierarchy();
		final PartitionCursor left=new PartitionCursor(), right=new PartitionCursor();
		int active=n;
		capture(result, groups, labels, targets, active);
		while(active>minimum){
			final String[] sorted=texts(groups);
			double best=-1;
			int bestA=-1, bestB=-1;
			for(int a=0; a<n; a++){
				if(groups[a]==null){continue;}
				for(int b=a+1; b<n; b++){
					if(groups[b]==null){continue;}
					final double value=scores[a*n+b];
					if(value<best){continue;}
					if(value>best || bestA<0 || before(groups, labels, separator, mergedCache, sorted, a, b, bestA, bestB, left, right)){
						best=value; bestA=a; bestB=b;
					}
				}
			}
			assert(bestA>=0 && bestB>bestA) : "At least two active groups supply a merge candidate";
			final Group a=groups[bestA], b=groups[bestB];
			final String merged=mergedText(groups, labels, separator, mergedCache, bestA, bestB);
			for(int c=0; c<n; c++){
				if(c!=bestA && c!=bestB && groups[c]!=null){
					scores[bestA*n+c]=scores[c*n+bestA]=Math.min(scores[bestA*n+c], scores[bestB*n+c]);
				}
				mergedCache[bestA*n+c]=mergedCache[c*n+bestA]=null;
				mergedCache[bestB*n+c]=mergedCache[c*n+bestB]=null;
			}
			groups[bestA]=new Group(union(a.members, b.members), merged); groups[bestB]=null;
			active--;
			result.merges.add(new Merge(active, a.text, b.text, merged, best));
			capture(result, groups, labels, targets, active);
		}
		return result;
	}

	private static boolean before(Group[] groups, String[] labels, String separator, String[] cache, String[] sorted,
			int a, int b, int bestA, int bestB, PartitionCursor left, PartitionCursor right){
		left.reset(sorted, groups[a].text, groups[b].text, mergedText(groups, labels, separator, cache, a, b));
		right.reset(sorted, groups[bestA].text, groups[bestB].text, mergedText(groups, labels, separator, cache, bestA, bestB));
		int x, y;
		do{x=left.next(); y=right.next(); if(x!=y){return x<y;}}while(x>=0);
		return false;
	}

	private static String mergedText(Group[] groups, String[] labels, String separator, String[] cache, int a, int b){
		final int slot=a*groups.length+b;
		if(cache[slot]==null){
			final StringBuilder text=new StringBuilder();
			for(int id : union(groups[a].members, groups[b].members)){
				if(text.length()>0){text.append(separator);} text.append(labels[id]);
			}
			cache[slot]=text.toString();
		}
		return cache[slot];
	}

	private static int[] union(int[] a, int[] b){
		final int[] out=new int[a.length+b.length];
		int i=0, j=0;
		for(int k=0; k<out.length; k++){
			out[k]=(j==b.length || (i<a.length && a[i]<b[j])) ? a[i++] : b[j++];
			assert(k==0 || out[k]>out[k-1]) : "Active groups are disjoint, sorted partitions of the dictionary";
		}
		return out;
	}

	private static String[] texts(Group[] groups){
		int count=0; for(Group group : groups){if(group!=null){count++;}}
		final String[] out=new String[count]; int pos=0;
		for(Group group : groups){if(group!=null){out[pos++]=group.text;}}
		Arrays.sort(out);
		return out;
	}

	private static void capture(Hierarchy result, Group[] groups, String[] labels, int[] targets, int active){
		boolean requested=false;
		for(int target : targets){if(target==active){requested=true;}}
		if(!requested){return;}
		final String[] text=texts(groups);
		final int[] map=new int[labels.length]; Arrays.fill(map, -1);
		for(Group group : groups){
			if(group!=null){
				final int id=Arrays.binarySearch(text, group.text);
				for(int symbol : group.members){require(map[symbol]==-1, "A symbol cannot belong to two reduced classes"); map[symbol]=id;}
			}
		}
		for(int id : map){require(id>=0 && id<active, "Every original symbol must survive in exactly one class");}
		result.cuts.add(new Cut(active, map, text));
	}

	private static final class Group{
		Group(int[] members_, String text_){members=members_; text=text_;}
		final int[] members;
		final String text;
	}

	/** Allocation-free view of a hypothetical sorted, slash-joined partition. */
	private static final class PartitionCursor{
		void reset(String[] sorted_, String a_, String b_, String merged_){
			sorted=sorted_; a=a_; b=b_; merged=merged_; current=null;
			index=position=emitted=0; inserted=separator=false;
		}
		int next(){
			if(separator){separator=false; return '/';}
			if(current==null){current=group(); position=0; if(current==null){return -1;}}
			final int value=current.charAt(position++);
			if(position==current.length()){current=null; emitted++; separator=emitted<sorted.length-1;}
			return value;
		}
		private String group(){
			while(index<sorted.length && (sorted[index]==a || sorted[index]==b)){index++;}
			if(!inserted && (index==sorted.length || merged.compareTo(sorted[index])<0)){inserted=true; return merged;}
			return index<sorted.length ? sorted[index++] : null;
		}
		String[] sorted;
		String a, b, merged, current;
		int index, position, emitted;
		boolean inserted, separator;
	}

	public static final class Hierarchy{
		public final ArrayList<Merge> merges=new ArrayList<Merge>();
		public final ArrayList<Cut> cuts=new ArrayList<Cut>();
	}
	public static final class Merge{
		Merge(int symbols_, String a_, String b_, String joined_, double odds_){symbols=symbols_; a=a_; b=b_; joined=joined_; odds=odds_;}
		public final int symbols;
		public final String a, b, joined;
		public final double odds;
	}
	public static final class Cut{
		Cut(int symbols_, int[] map_, String[] groups_){symbols=symbols_; map=map_; groups=groups_;}
		public final int symbols;
		public final int[] map;
		public final String[] groups;
	}

	private static long sum(long[] values){
		long total=0;
		for(long value : values){require(value>=0, "Observation counters must not be negative or overflowed"); total=Math.addExact(total, value);}
		return total;
	}
	private static void require(boolean condition, String why){if(!condition){throw new IllegalArgumentException(why);}}
}
