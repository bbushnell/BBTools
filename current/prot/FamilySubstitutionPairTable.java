package prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicInteger;

import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Empirical amino-acid substitution tables from right-pinned family-member tails.
 * The input is a pinned FamilyCalibration.MemberInputs preparation. Each member is
 * aligned to its accepted family consensus with the suffix geometry used by the
 * repaired CallGenes miss path: the consensus is consumed globally, the member
 * consumes a nonempty suffix, and only aligned standard-residue columns are counted.
 *
 * @author Yelan
 */
public final class FamilySubstitutionPairTable{

	private FamilySubstitutionPairTable(){}

	public static void main(String[] args){
		try{
			final PreParser pp=new PreParser(args, FamilySubstitutionPairTable.class, false);
			try{
				if(pp.args.length==1 && pp.args[0].equalsIgnoreCase("test")){test();}
				else{run(pp.args);}
			}finally{Shared.closeStream(pp.outstream);}
		}catch(Throwable failure){failure.printStackTrace(); System.exit(1);}
	}

	public static void run(String[] args) throws Exception{
		final HashMap<String,String> options=parse(args);
		final String population=population(options);
		final Path out=Paths.get(required(options, "out")).toAbsolutePath().normalize();
		final boolean overwrite=options.containsKey("overwrite") && Parse.parseBoolean(options.get("overwrite"));
		require((overwrite || !Files.exists(out)) && Files.isDirectory(out.getParent()),
			"Output directory must be fresh unless overwrite=t and under an existing parent: "+out);
		if(!Files.exists(out)){Files.createDirectory(out);}

		final int threads=options.containsKey("t") ? integer(options, "t") : Shared.threads();
		final int expectedFamilies=integer(options, "expectedfamilies");
		final long expectedMembers=number(options, "expectedmembers");
		final int buckets=integer(options, "buckets");
		final int batches=options.containsKey("batches") ? integer(options, "batches") : integer(options, "expectedbatches");
		final int maxReference=integer(options, "maxreferenceaa");
		final long maxMembers=number(options, "maxbucketmembers"), maxBytes=number(options, "maxbucketbytes");
		final IdentityBands bands=options.containsKey("bands") ? IdentityBands.parse(options.get("bands")) : IdentityBands.legacy();
		final FamilyCalibration.VerificationMode verification=options.containsKey("verification") ?
			FamilyCalibration.VerificationMode.valueOf(options.get("verification").toUpperCase(Locale.ROOT)) :
			FamilyCalibration.VerificationMode.SHA80;
		require(threads>0 && expectedFamilies>0 && expectedMembers>=0 && buckets>0 && batches>0
			&& maxReference>0 && maxMembers>0 && maxBytes>0, "Invalid thread/count/cap dimensions");

		Shared.AMINO_IN=true; Shared.setThreads(threads); ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=false;
		final FamilyCalibration.MemberInputs inputs=new FamilyCalibration.MemberInputs(
			input(options, "families"), required(options, "familiessha80"),
			input(options, "chunks"), required(options, "chunkssha80"),
			input(options, "consensus"), required(options, "consensussha80"),
			expectedFamilies, expectedMembers, buckets, batches, maxReference, MemberLengthRange.from(options), verification);
		final Counts total=new Counts(bands.size());
		final long start=System.nanoTime();
		for(int bucket=0; bucket<inputs.buckets; bucket++){
			final FamilyCalibration.MemberInput[] families=inputs.loadBucket(bucket, maxMembers, maxBytes);
			final Counts[] partial=process(families, threads, bands);
			for(Counts counts : partial){total.add(counts);}
		}
		inputs.verifyComplete();
		writeOutputs(out, total, args, expectedFamilies, expectedMembers, start, System.nanoTime(), population, bands, verification);
		System.out.println("FAMILY_SUBSTITUTION_PAIR_TABLE_PASS families="+total.families+
			" members="+total.members+" aligned_pairs="+total.alignedPairs());
	}

	private static Counts[] process(FamilyCalibration.MemberInput[] families, int threads, IdentityBands bands) throws InterruptedException{
		final int workers=Math.min(threads, families.length);
		final Counts[] out=new Counts[workers];
		final AtomicInteger next=new AtomicInteger();
		final Worker[] array=new Worker[workers];
		for(int i=0; i<workers; i++){
			out[i]=new Counts(bands.size());
			array[i]=new Worker(i, families, next, out[i], bands);
		}
		for(Worker worker : array){worker.start();}
		for(Worker worker : array){worker.join();}
		for(Worker worker : array){if(worker.failure!=null){throw new IllegalStateException("Substitution worker failed: "+worker.getName(), worker.failure);}}
		return out;
	}

	private static final class Worker extends Thread{
		Worker(int id, FamilyCalibration.MemberInput[] families_, AtomicInteger next_, Counts counts_, IdentityBands bands_){
			super("family-substitution-"+id);
			families=families_; next=next_; counts=counts_; bands=bands_;
		}
		@Override public void run(){
			try{
				for(int index=next.getAndIncrement(); index<families.length; index=next.getAndIncrement()){
					process(families[index]);
					families[index]=null;
				}
			}catch(Throwable caught){failure=caught; next.set(families.length);}
		}
		private void process(FamilyCalibration.MemberInput family){
			final ProteinSequence consensus=new ProteinSequence(family.id, family.consensus);
			final Counts local=new Counts(bands.size());
			counts.families++;
			counts.members+=family.members.length;
			for(byte[] raw : family.members){
				final ProteinSequence member=buffers.encode(raw);
				final AAAlignment alignment=StopAnchoredAminoLinear.align(member.enc, consensus.enc,
					StopAnchoredAminoLinear.Mode.QUERY_SUFFIX, true);
				final int mask=bands.mask(alignment.pident());
				walk(member.enc, consensus.enc, alignment, mask, counts);
				walk(member.enc, consensus.enc, alignment, mask, local);
				for(int band=0; band<bands.size(); band++){
					if((mask&(1<<band))!=0){
						counts.alignments[band]++;
						local.alignments[band]++;
					}
				}
			}
			counts.addFamilyEqual(local);
		}
		private final FamilyCalibration.MemberInput[] families;
		private final AtomicInteger next;
		private final Counts counts;
		private final IdentityBands bands;
		private final FamilyCalibration.Buffers buffers=new FamilyCalibration.Buffers();
		Throwable failure;
	}

	private static void walk(byte[] query, byte[] target, AAAlignment alignment, int mask, Counts counts){
		if(mask==0){return;}
		int q=alignment.qStart, t=alignment.tStart;
		for(byte op : alignment.match){
			if(op=='m'){
				final int qa=query[q++], ta=target[t++];
				if(Blosum62.isStandard((byte)qa) && Blosum62.isStandard((byte)ta)){
					for(int band=0; band<counts.alignments.length; band++){
						if((mask&(1<<band))!=0){
							counts.ordered[band][qa][ta]++;
							counts.query[band][qa]++;
							counts.consensus[band][ta]++;
							counts.combined[band][qa]++;
							counts.combined[band][ta]++;
							counts.symmetric[band][Math.min(qa, ta)][Math.max(qa, ta)]++;
						}
					}
				}
			}else if(op=='I'){q++;}
			else if(op=='D'){t++;}
			else{throw new AssertionError("Unknown alignment op: "+(char)op);}
		}
	}

	private static void writeOutputs(Path out, Counts counts, String[] args, int expectedFamilies,
			long expectedMembers, long start, long stop, String population, IdentityBands bands,
			FamilyCalibration.VerificationMode verification) throws Exception{
		writeOrdered(out.resolve("ordered_counts.tsv"), counts, bands);
		final Table[] tables=writeSymmetric(out.resolve("symmetric_tables.tsv"), counts, bands);
		writeGroups(out.resolve("merge_groups.tsv"), tables);
		write(out.resolve("summary.tsv"), "families\tmembers\talignments\taligned_pairs\tseconds\n"+
			counts.families+"\t"+counts.members+"\t"+counts.alignments()+"\t"+counts.alignedPairs()+"\t"+((stop-start)*1e-9)+"\n");
		write(out.resolve("recipe.txt"), "schema=family_substitution_pair_table_v1\n"+
			"population="+population+"\n"+
			"verification="+verification+"\n"+
			"alignment=StopAnchoredAminoLinear_QUERY_SUFFIX_right_pinned_left_clipped\n"+
			"bands="+bands.recipe()+"\n"+
			"expectedfamilies="+expectedFamilies+"\nexpectedmembers="+expectedMembers+"\n"+
			String.join("\n", args)+"\n");
		write(out.resolve("PASS"), "FAMILY_SUBSTITUTION_PAIR_TABLE_PASS\n");
	}

	private static void writeOrdered(Path path, Counts counts, IdentityBands bands){
		final ByteStreamWriter writer=writer(path);
		final ByteBuilder line=new ByteBuilder();
		try{
			writer.println("band\tquery_aa\tconsensus_aa\tcount\tquery_monomer\tconsensus_monomer\talignments\taligned_pairs");
			for(int band=0; band<bands.size(); band++){
				final long aligned=counts.alignedPairs(band);
				for(int q=0; q<AA_COUNT; q++){
					for(int t=0; t<AA_COUNT; t++){
						line.clear().append(bands.label(band)).tab().append((char)AminoAcid.numberToAcid[q]).tab()
							.append((char)AminoAcid.numberToAcid[t]).tab().append(counts.ordered[band][q][t]).tab()
							.append(counts.query[band][q]).tab().append(counts.consensus[band][t]).tab()
							.append(counts.alignments[band]).tab().append(aligned);
						writer.println(line);
					}
				}
			}
		}finally{require(!writer.poisonAndWait(), "Ordered count output failed");}
	}

	private static Table[] writeSymmetric(Path path, Counts counts, IdentityBands bands){
		final Table[] tables=new Table[bands.size()*2];
		final ByteStreamWriter writer=writer(path);
		final ByteBuilder line=new ByteBuilder();
		try{
			writer.println("weighting\tband\taa_i\taa_j\tcount\tobserved_frequency\texpected_frequency\todds\tfamilies\talignments\taligned_pairs");
			int table=0;
			for(int band=0; band<bands.size(); band++){
				tables[table]=Table.raw(bands.label(band), counts.symmetric[band], counts.combined[band]);
				final long alignedPairs=counts.alignedPairs(band);
				writeTable(writer, line, "raw", bands.label(band), tables[table++], 0, counts.alignments[band], alignedPairs);
				tables[table]=Table.familyEqual(bands.label(band), counts.familySymmetric[band], counts.familyCombined[band],
					counts.familyBand[band]);
				writeTable(writer, line, "family_equal", bands.label(band), tables[table++],
					counts.familyBand[band], counts.alignments[band], alignedPairs);
			}
		}finally{require(!writer.poisonAndWait(), "Symmetric table output failed");}
		return tables;
	}

	private static void writeTable(ByteStreamWriter writer, ByteBuilder line, String weighting, String band,
			Table table, long families, long alignments, long alignedPairs){
		for(int i=0; i<AA_COUNT; i++){
			for(int j=i; j<AA_COUNT; j++){
				line.clear().append(weighting).tab().append(band).tab().append((char)AminoAcid.numberToAcid[i]).tab()
					.append((char)AminoAcid.numberToAcid[j]).tab().appendSlow(table.count[i][j]).tab()
					.append(Double.toString(table.observed[i][j])).tab().append(Double.toString(table.expected[i][j])).tab()
					.append(Double.toString(table.odds[i][j])).tab().append(families).tab().append(alignments).tab()
					.append(alignedPairs);
				writer.println(line);
			}
		}
	}

	private static void writeGroups(Path path, Table[] tables){
		final ByteStreamWriter writer=writer(path);
		final ByteBuilder line=new ByteBuilder();
		try{
			writer.println("weighting\tband\trule\tsymbols\tgroups");
			for(int i=0; i<tables.length/2; i++){
				writeGroups(writer, line, "raw", tables[i*2]);
				writeGroups(writer, line, "family_equal", tables[i*2+1]);
			}
			for(String[] row : REFERENCE_GROUPS){
				line.clear().append(row[0]).tab().append("reference").tab().append("published").tab().append(row[1]).tab().append(row[2]);
				writer.println(line);
			}
		}finally{require(!writer.poisonAndWait(), "Merge group output failed");}
	}

	private static void writeGroups(ByteStreamWriter writer, ByteBuilder line, String weighting, Table table){
		for(String rule : new String[]{"weakest_link_odds", "seed_objective"}){
			for(int symbols : new int[]{12, 10, 8}){
				line.clear().append(weighting).tab().append(table.band).tab().append(rule).tab().append(symbols).tab()
					.append(reduce(table, rule, symbols));
				writer.println(line);
			}
		}
	}

	private static String reduce(Table table, String rule, int target){
		final ArrayList<Integer> groups=new ArrayList<Integer>(AA_COUNT);
		for(int i=0; i<AA_COUNT; i++){groups.add(1<<i);}
		while(groups.size()>target){
			double best=Double.NEGATIVE_INFINITY;
			int bestA=-1, bestB=-1;
			String bestText=null;
			for(int a=0; a<groups.size(); a++){
				for(int b=a+1; b<groups.size(); b++){
					final double score=rule.equals("weakest_link_odds") ?
						minOdds(table, groups.get(a), groups.get(b)) : seedObjective(table, groups.get(a), groups.get(b));
					final String text=format(merged(groups, a, b));
					if(score>best || (score==best && (bestText==null || text.compareTo(bestText)<0))){
						best=score; bestA=a; bestB=b; bestText=text;
					}
				}
			}
			groups.set(bestA, groups.get(bestA)|groups.get(bestB));
			groups.remove(bestB);
		}
		return format(groups);
	}

	private static double minOdds(Table table, int a, int b){
		double min=Double.POSITIVE_INFINITY;
		for(int i=0; i<AA_COUNT; i++){
			if((a&(1<<i))==0){continue;}
			for(int j=0; j<AA_COUNT; j++){
				if((b&(1<<j))!=0){min=Math.min(min, table.odds[Math.min(i, j)][Math.max(i, j)]);}
			}
		}
		return min;
	}

	private static double seedObjective(Table table, int a, int b){
		double observed=0, expected=0;
		for(int i=0; i<AA_COUNT; i++){
			if((a&(1<<i))==0){continue;}
			for(int j=0; j<AA_COUNT; j++){
				if((b&(1<<j))!=0){
					observed+=table.observed[Math.min(i, j)][Math.max(i, j)];
					expected+=table.expected[Math.min(i, j)][Math.max(i, j)];
				}
			}
		}
		return expected>0 ? observed/expected : -1;
	}

	private static ArrayList<Integer> merged(ArrayList<Integer> groups, int a, int b){
		final ArrayList<Integer> copy=new ArrayList<Integer>(groups);
		copy.set(a, copy.get(a)|copy.get(b));
		copy.remove(b);
		return copy;
	}

	private static String format(ArrayList<Integer> groups){
		final ArrayList<String> text=new ArrayList<String>(groups.size());
		for(int mask : groups){
			final StringBuilder sb=new StringBuilder();
			for(int i=0; i<AA_COUNT; i++){
				if((mask&(1<<i))!=0){sb.append((char)AminoAcid.numberToAcid[i]);}
			}
			text.add(sb.toString());
		}
		java.util.Collections.sort(text);
		final StringBuilder sb=new StringBuilder();
		for(String s : text){
			if(sb.length()>0){sb.append('/');}
			sb.append(s);
		}
		return sb.toString();
	}

	private static void test() throws Exception{
		final IdentityBands legacy=IdentityBands.legacy();
		final Counts counts=new Counts(legacy.size()), family=new Counts(legacy.size());
		final ProteinSequence query=new ProteinSequence("q", "ACDEFGHIKLMNPQRSTVWY");
		final ProteinSequence target=new ProteinSequence("t", "AYDEFGHIKLMNPQRSTVWY");
		final AAAlignment alignment=StopAnchoredAminoLinear.align(query.enc, target.enc,
			StopAnchoredAminoLinear.Mode.QUERY_SUFFIX, true);
		final int band=legacy.exclusiveBand(alignment.pident());
		walk(query.enc, target.enc, alignment, 1<<band, counts);
		walk(query.enc, target.enc, alignment, 1<<band, family);
		counts.alignments[band]++; family.alignments[band]++;
		counts.families=1; counts.members=1; counts.addFamilyEqual(family);
		require(counts.ordered[2][AminoAcid.acidToNumber['A']][AminoAcid.acidToNumber['A']]==1, "Missing A/A test pair");
		require(counts.ordered[2][AminoAcid.acidToNumber['C']][AminoAcid.acidToNumber['Y']]==1, "Missing ordered C/Y test pair");
		require(counts.symmetric[2][AminoAcid.acidToNumber['C']][AminoAcid.acidToNumber['Y']]==1, "Missing symmetric C/Y test pair");
		require(counts.alignedPairs()==20 && counts.alignments()==1, "Bad test totals");
		final Table table=Table.raw(">=60", counts.symmetric[2], counts.combined[2]);
		require(reduce(table, "weakest_link_odds", 12).indexOf('/')>=0, "Merge ladder failed");
		testPopulation(counts);
		final IdentityBands custom=IdentityBands.parse("0to40,25to45");
		require(custom.mask(30)==3 && custom.mask(50)==0 && custom.mask(25)==3 && custom.mask(45)==0,
			"Overlapping custom identity bands must be lower-inclusive and upper-exclusive");
		System.out.println("FAMILY_SUBSTITUTION_PAIR_TABLE_TEST_PASS");
	}

	/** Checks emitted provenance and byte-identical numerical outputs for an explicit population. */
	private static void testPopulation(Counts counts) throws Exception{
		final String baseline=population(parse(new String[]{"batches=1"}));
		require(baseline.equals("all_existing_K5_minhits6_winner_assigned_members"),
			"Omitting population must preserve the accepted baseline provenance bytes");
		final String label="fixed_cap150_K4p1_minhits5_recruited_members";
		final String[] args={"batches=1", "population="+label};
		final Path root=Files.createTempDirectory("family-substitution-provenance-");
		final Path original=Files.createDirectory(root.resolve("default"));
		final Path explicit=Files.createDirectory(root.resolve("explicit"));
		writeOutputs(original, counts, new String[]{"batches=1"}, 1, 1, 0, 0, baseline, IdentityBands.legacy(),
			FamilyCalibration.VerificationMode.SHA80);
		writeOutputs(explicit, counts, args, 1, 1, 0, 0, population(parse(args)), IdentityBands.legacy(),
			FamilyCalibration.VerificationMode.SHA80);
		final String oldText=new String(Files.readAllBytes(original.resolve("recipe.txt")), StandardCharsets.UTF_8);
		final String newText=new String(Files.readAllBytes(explicit.resolve("recipe.txt")), StandardCharsets.UTF_8);
		require(oldText.startsWith("schema=family_substitution_pair_table_v1\npopulation="+baseline+"\n"),
			"Default output must identify the original K5/minhits6 member population");
		require(newText.startsWith("schema=family_substitution_pair_table_v1\npopulation="+label+"\n"),
			"Explicit population must reach recipe.txt verbatim rather than mislabel recruited members");
		for(String file : new String[]{"ordered_counts.tsv", "symmetric_tables.tsv", "merge_groups.tsv", "summary.tsv", "PASS"}){
			require(Arrays.equals(Files.readAllBytes(original.resolve(file)), Files.readAllBytes(explicit.resolve(file))),
				"Population metadata must not alter numerical output: "+file);
		}
		for(String invalid : new String[]{"", "bad\tlabel", "bad\nlabel", "bad\rlabel"}){
			boolean rejected=false;
			try{population(parse(new String[]{"batches=1", "population="+invalid}));}
			catch(IllegalArgumentException expected){rejected=true;}
			require(rejected, "Population labels must not inject empty values or extra recipe fields");
		}
	}

	/** A population is provenance, never a scoring parameter; preserve the historical default. */
	private static String population(HashMap<String,String> options){
		final String label=options.containsKey("population") ? required(options, "population") : DEFAULT_POPULATION;
		for(int i=0; i<label.length(); i++){
			require(!Character.isISOControl(label.charAt(i)), "Population label must remain one printable recipe field");
		}
		return label;
	}

	private static HashMap<String,String> parse(String[] args){
		final HashMap<String,String> options=new HashMap<String,String>();
		for(String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected flag=value: "+arg);
			String key=pair[0].toLowerCase(Locale.ROOT), value=pair[1];
			if(key.equals("threads")){key="t";}
			else if(key.equals("ow")){key="overwrite";}
			require(OPTIONS.contains(key) && options.put(key, value)==null, "Unknown or duplicate option: "+key);
		}
		require(options.containsKey("batches") || options.containsKey("expectedbatches"), "Require batches= or expectedbatches=");
		return options;
	}

	private static Path input(HashMap<String,String> o, String key){
		return Paths.get(required(o, key)).toAbsolutePath().normalize();
	}
	private static String required(HashMap<String,String> o, String key){
		final String value=o.get(key);
		require(value!=null && !value.isEmpty(), "Missing option: "+key);
		return value;
	}
	private static long number(HashMap<String,String> o, String key){
		final long value=Long.parseLong(required(o, key));
		require(value>=0, "Negative option: "+key);
		return value;
	}
	private static int integer(HashMap<String,String> o, String key){return Math.toIntExact(number(o, key));}
	private static void write(Path path, String text) throws Exception{Files.write(path, text.getBytes(StandardCharsets.UTF_8));}
	private static ByteStreamWriter writer(Path path){
		final ByteStreamWriter writer=new ByteStreamWriter(path.toString(), true, false, false);
		writer.start();
		return writer;
	}
	private static void require(boolean condition, String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	private static final class Counts{
		Counts(int bands){
			alignments=new long[bands]; familyBand=new long[bands];
			ordered=new long[bands][AA_COUNT][AA_COUNT];
			symmetric=new long[bands][AA_COUNT][AA_COUNT];
			query=new long[bands][AA_COUNT]; consensus=new long[bands][AA_COUNT]; combined=new long[bands][AA_COUNT];
			familySymmetric=new double[bands][AA_COUNT][AA_COUNT];
			familyCombined=new double[bands][AA_COUNT];
		}
		void add(Counts b){
			families+=b.families; members+=b.members;
			for(int band=0; band<alignments.length; band++){
				alignments[band]+=b.alignments[band]; familyBand[band]+=b.familyBand[band];
				for(int i=0; i<AA_COUNT; i++){
					query[band][i]+=b.query[band][i]; consensus[band][i]+=b.consensus[band][i]; combined[band][i]+=b.combined[band][i];
					familyCombined[band][i]+=b.familyCombined[band][i];
					for(int j=0; j<AA_COUNT; j++){ordered[band][i][j]+=b.ordered[band][i][j];}
					for(int j=i; j<AA_COUNT; j++){
						symmetric[band][i][j]+=b.symmetric[band][i][j];
						familySymmetric[band][i][j]+=b.familySymmetric[band][i][j];
					}
				}
			}
		}
		void addFamilyEqual(Counts family){
			for(int band=0; band<alignments.length; band++){
				final long total=family.alignedPairs(band);
				if(total<=0){continue;}
				familyBand[band]++;
				final double pairDen=1.0/total, monoDen=0.5/total;
				for(int i=0; i<AA_COUNT; i++){
					familyCombined[band][i]+=family.combined[band][i]*monoDen;
					for(int j=i; j<AA_COUNT; j++){familySymmetric[band][i][j]+=family.symmetric[band][i][j]*pairDen;}
				}
			}
		}
		long alignments(){long sum=0; for(long v : alignments){sum+=v;} return sum;}
		long alignedPairs(){long sum=0; for(int b=0; b<alignments.length; b++){sum+=alignedPairs(b);} return sum;}
		long alignedPairs(int band){long sum=0; for(int i=0; i<AA_COUNT; i++){for(int j=i; j<AA_COUNT; j++){sum+=symmetric[band][i][j];}} return sum;}

		long families, members;
		final long[] alignments, familyBand;
		final long[][][] ordered, symmetric;
		final long[][] query, consensus, combined;
		final double[][][] familySymmetric;
		final double[][] familyCombined;
	}

	private static final class Table{
		static Table raw(String band, long[][] pairs, long[] monomers){
			final Table table=new Table(band);
			for(int i=0; i<AA_COUNT; i++){for(int j=i; j<AA_COUNT; j++){table.total+=pairs[i][j];}}
			if(table.total>0){
				for(int i=0; i<AA_COUNT; i++){table.monomer[i]=monomers[i]/(2.0*table.total);}
				for(int i=0; i<AA_COUNT; i++){
					for(int j=i; j<AA_COUNT; j++){
						table.count[i][j]=pairs[i][j];
						table.observed[i][j]=pairs[i][j]/table.total;
					}
				}
			}
			table.finish();
			return table;
		}
		static Table familyEqual(String band, double[][] pairs, double[] monomers, long families){
			final Table table=new Table(band);
			if(families>0){
				table.total=families;
				for(int i=0; i<AA_COUNT; i++){table.monomer[i]=monomers[i]/families;}
				for(int i=0; i<AA_COUNT; i++){
					for(int j=i; j<AA_COUNT; j++){
						table.count[i][j]=pairs[i][j];
						table.observed[i][j]=pairs[i][j]/families;
					}
				}
			}
			table.finish();
			return table;
		}
		private Table(String band_){band=band_;}
		private void finish(){
			for(int i=0; i<AA_COUNT; i++){
				for(int j=i; j<AA_COUNT; j++){
					expected[i][j]=(i==j ? monomer[i]*monomer[i] : 2*monomer[i]*monomer[j]);
					odds[i][j]=expected[i][j]>0 ? observed[i][j]/expected[i][j] : 0;
				}
			}
		}
		final String band;
		double total;
		final double[] monomer=new double[AA_COUNT];
		final double[][] count=new double[AA_COUNT][AA_COUNT];
		final double[][] observed=new double[AA_COUNT][AA_COUNT], expected=new double[AA_COUNT][AA_COUNT], odds=new double[AA_COUNT][AA_COUNT];
	}

	private static final int AA_COUNT=20;
	private static final String DEFAULT_POPULATION="all_existing_K5_minhits6_winner_assigned_members";
	private static final ArrayList<String> OPTIONS=new ArrayList<String>(Arrays.asList(
		"families", "familiessha80", "chunks", "chunkssha80", "consensus", "consensussha80",
		"expectedfamilies", "expectedmembers", "expectedbatches", "batches", "buckets",
		"maxreferenceaa", "maxbucketmembers", "maxbucketbytes", "out", "t", "overwrite",
		"minmemberaa", "maxmemberaa", "population", "bands", "verification"));
	private static final String[][] REFERENCE_GROUPS={
		{"diamond11", "11", "KREDQN/C/G/H/ILV/M/F/Y/W/P/STA"},
		{"steinegger12_mmseqs", "12", "AST/C/DN/EQ/FY/G/H/IV/KR/LM/P/W"},
	};
}
