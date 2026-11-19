package benchmark;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.LinkOption;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.attribute.BasicFileAttributes;
import java.util.HashMap;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import parse.Parse;
import parse.LineParser1;
import parse.PreParser;
import prot.DigestSuffix;
import prot.FamilyReferenceCodec;
import prot.FamilyReference.AlignmentCalibration;
import prot.FamilyGateCodec;
import prot.FamilyTailCascade;
import prot.MemberLengthRange;
import prot.FamilyCalibration.MemberInput;
import prot.FamilyCalibration.MemberInputs;
import prot.FamilyCalibration.VerificationMode;
import shared.Shared;
import structures.ByteBuilder;

/** Counts natural member lengths through the canonical pinned preparation reader.
 * No sequence classification, alignment, family assignment or consensus update occurs.
 * @author Fischl
 */
public final class FamilyMemberLengthCensus{

	public static void main(final String[] args) throws Exception{
		final PreParser parser=new PreParser(args, FamilyMemberLengthCensus.class, false);
		try{
			final HashMap<String,String> options=new HashMap<String,String>();
			for(String arg:parser.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				require(pair.length==2 && (pair[0].equals("in") || pair[0].equals("insha80") || pair[0].equals("out")
					|| pair[0].equals("histogram") || pair[0].equals("members") || pair[0].equals("expected")
					|| pair[0].equals("reference") || pair[0].equals("marginscale") || pair[0].equals("verification")),
					"Expected preparation in/insha80 or histogram/members/expected, plus out=new_directory");
				require(options.put(pair[0], pair[1])==null, "Duplicate census option: "+pair[0]);
			}
			if(options.containsKey("reference")){
				require(options.size()==3+(options.containsKey("insha80") ? 1 : 0)+(options.containsKey("verification") ? 1 : 0),
					"Threshold mode needs reference, marginscale, out and optional insha80/verification");
				thresholds(Paths.get(required(options, "reference")), options.getOrDefault("insha80", "unmeasured"),
					Double.parseDouble(required(options, "marginscale")), Paths.get(required(options, "out")),
					verification(options, "threshold"));
			}else if(options.containsKey("histogram")){
				require(options.size()==4, "Coverage mode needs exactly histogram, members, expected and out");
				coverage(Paths.get(required(options, "histogram")), Paths.get(required(options, "members")),
					Long.parseLong(required(options, "expected")), Paths.get(required(options, "out")));
			}else{
				require(options.size()==2+(options.containsKey("insha80") ? 1 : 0)+(options.containsKey("verification") ? 1 : 0),
					"Member mode needs in, out and optional insha80/verification");
				count(Paths.get(required(options, "in")), options.getOrDefault("insha80", "unmeasured"),
					Paths.get(required(options, "out")), verification(options, "member"));
			}
		}finally{Shared.closeStream(parser.outstream);}
	}

	/*--------------------------------------------------------------*/

	/** Reports measured saved gates, with explicit unavailable and empty-family counts.
	 * Header-bound alignment geometry is reported; payload lines are skipped, never recalibrated.
	 */
	private static void thresholds(final Path reference, final String pin, final double scale, final Path out,
			final VerificationMode verification) throws Exception{
		final BasicFileAttributes initial;
		if(verification==VerificationMode.STRUCTURAL){
			require(pin.equals("unmeasured") || pin.matches("[0-9a-f]{20}"), "Invalid structural reference identity: "+pin);
			initial=structuralState(reference);
		}else{
			require(pin.matches("[0-9a-f]{20}") && pin.equals(DigestSuffix.file(reference.toString())), "Reference fingerprint differs");
			initial=null;
		}
		require(Double.isFinite(scale) && scale>=0 && !Files.exists(out), "Threshold summary needs a valid margin and fresh output");
		final ByteFile input=ByteFile.makeByteFile1(reference.toString(), false);
		final LineParser1 parser=new LineParser1('\t');
		final byte[][] block=new byte[5][];
		final double[][] values;
		final FamilyReferenceCodec.Header header;
		int live=0, measured=0, unavailable=0, empty=0;
		try{
			for(int i=0; i<block.length; i++){block[i]=input.nextLine();}
			header=FamilyReferenceCodec.parseMetricBoundHeader(block, 0, parser);
			require(header.liveFamilies>0 && header.totalFamilies<=1000000, "Bounded nonempty reference roster required");
			// The canonical codec validated this exact schema's labels, including the optional state column.
			parser.set(block[1]);
			final int gateTerms=parser.terms()-1;
			require(gateTerms==17 || gateTerms==18, "Validated reference header must name the saved gate layout");
			values=new double[3][header.liveFamilies];
			final boolean[] seen=new boolean[header.totalFamilies];
			for(int row=0; row<header.totalFamilies; row++){
				for(int i=0; i<block.length; i++){
					block[i]=input.nextLine();
					require(block[i]!=null, "Every family needs all five saved reference lines");
				}
				parser.set(block[0]);
				require(parser.terms()==6 && parser.termEquals("#family", 0), "Expected native family metadata");
				final int rank=parser.parseInt(2);
				final long members=parser.parseLong(3);
				require(rank>=0 && rank<seen.length && !seen[rank] && members>=0, "Each physical rank must appear once");
				seen[rank]=true;
				if(members==0){
					require(parser.parseLong(4)==0 && parser.parseLong(5)==0,
						"Empty family must have zero residue and word counts");
					for(int i=1; i<block.length; i++){
						require(block[i].length==2 && block[i][0]=='N' && block[i][1]=='A', "Empty family payloads must be NA");
					}
					empty++;
					continue;
				}
				parser.set(block[1]);
				require(parser.terms()==gateTerms && live<header.liveFamilies,
					"Saved gate fields must match the header and bounded live roster");
				live++;
				final boolean enabled=gateTerms==17 || parser.termEquals(AlignmentCalibration.ENABLED.label, 7);
				require(enabled || parser.termEquals(AlignmentCalibration.UNAVAILABLE.label, 7),
					"Unknown alignment calibration availability");
				final float lower=FamilyGateCodec.readLower(parser, 1), upper=FamilyGateCodec.readUpper(parser, 4);
				if(!enabled){
					require(lower==0 && upper==0 && FamilyGateCodec.readLower(parser, 2)==0 &&
						FamilyGateCodec.readFinalBlosum(parser, 5)==0,
						"Unavailable alignment gates must contain inert zeros, never measured thresholds");
					unavailable++;
					continue;
				}
				require(lower>=0 && lower<=100 && upper>=lower, "Saved identity gates must preserve member lower");
				values[0][measured]=lower; values[1][measured]=upper;
				values[2][measured]=FamilyTailCascade.evaluationGate(lower, upper, 30f, scale);
				measured++;
			}
			require(live==header.liveFamilies && empty==header.emptyFamilies && measured+unavailable==live &&
				input.nextLine()==null, "Threshold counts must partition the complete declared roster");
		}finally{require(!input.close(), "Reference read failed");}
		if(initial==null){require(pin.equals(DigestSuffix.file(reference.toString())), "Reference changed during summary");}
		else{verifyStable(reference, initial);}
		final String[] names={"member_identity_lower", "family_identity_upper", "effective_identity_accept_floor30"};
		final String status=unavailable==0 ? "ALL_LIVE_FAMILIES_MEASURED" : "UNEXPECTED_UNAVAILABLE";
		final ByteBuilder table=new ByteBuilder("metric\ttotal_families\tmeasured_families\tunavailable_families\tempty_families"
			+"\tavailability\tmargin_scale\tmin\tmedian\tp90\tmax\talignment_mode\tminimum_span\tminimum_paired\tmatrix_name\n");
		for(int i=0; i<values.length; i++){
			java.util.Arrays.sort(values[i], 0, measured);
			table.append(names[i]).tab().append(header.totalFamilies).tab().append(measured).tab()
				.append(unavailable).tab().append(empty).tab().append(status).tab().append(scale, 6);
			for(int percent:new int[]{0,50,90,100}){
				table.tab();
				if(measured==0){table.append("NA");}
				else{
					final int index=(int)Math.max(0, ((long)percent*measured+99)/100-1);
					table.append(values[i][index], 9);
				}
			}
			table.tab().append(header.alignmentMetric.name()).tab().append(header.alignmentMetric.minimumSpan)
				.tab().append(header.alignmentMetric.minimumPaired).tab().append(header.alignmentMetric.matrix.name).nl();
		}
		Files.createDirectory(out);
		final ByteStreamWriter writer=new ByteStreamWriter(out.resolve("thresholds.tsv").toString(), false, false, false);
		writer.start(); writer.print(table);
		require(!writer.poisonAndWait(), "Threshold table write failed");
		Files.write(out.resolve("PASS"), new ByteBuilder("FAMILY_IDENTITY_THRESHOLDS_PASS families=").append(live)
			.append(" measured=").append(measured).append(" unavailable=").append(unavailable).append(" empty=").append(empty)
			.append(" status=").append(status).append(" quantile=nearest_rank floor=30\n").toBytes());
		System.out.print(table.toString());
		System.err.println(unavailable==0 ? "Alignment threshold census: 0 unavailable families; all "+measured+
			" live families have measured thresholds." : "WARNING: "+unavailable+" of "+live+
			" live families have unavailable alignment thresholds; this is unexpected. Quantiles use only the "+measured+
			" measured families.");
	}
	/*----------------            Census            ----------------*/
	/*--------------------------------------------------------------*/

	/** Joins the native Reformat length histogram to the already verified member census. */
	private static void coverage(final Path histogram, final Path memberTable, final long expected, final Path out) throws Exception{
		require(expected>0 && !Files.exists(out), "Coverage requires a positive population count and fresh output");
		final long[] total=new long[4], assigned=new long[4];
		final boolean[] seen=new boolean[4], lengths=new boolean[201];
		final LineParser1 parser=new LineParser1('\t');
		final ByteFile hist=ByteFile.makeByteFile1(histogram.toString(), false);
		try{
			require(java.util.Arrays.equals(hist.nextLine(), "#Length\tCount".getBytes(StandardCharsets.US_ASCII)), "Expected native Reformat length histogram");
			for(byte[] line=hist.nextLine(); line!=null; line=hist.nextLine()){
				parser.set(line);
				require(parser.terms()==2, "Length histogram rows need length and count");
				final int length=parser.parseInt(0);
				final long count=parser.parseLong(1);
				require(count>=0 && length>=0 && length<=200 && !lengths[length], "Histogram lengths must be unique and at most200");
				lengths[length]=true;
				if(count>0){total[lengthBin(length)]=Math.addExact(total[lengthBin(length)], count);}
			}
		}finally{require(!hist.close(), "Length histogram read failed");}
		final ByteFile members=ByteFile.makeByteFile1(memberTable.toString(), false);
		long recordedMembers=-1;
		try{
			require(java.util.Arrays.equals(members.nextLine(), "metric\tcount".getBytes(StandardCharsets.US_ASCII)), "Expected member census table");
			for(byte[] line=members.nextLine(); line!=null; line=members.nextLine()){
				parser.set(line);
				require(parser.terms()==2, "Member census rows need metric and count");
				final long count=parser.parseLong(1);
				require(count>=0, "Member census cannot contain negative counts");
				if(parser.termEquals("members", 0)){
					require(recordedMembers<0, "Member total must occur exactly once");
					recordedMembers=count;
				}
				for(int i=0; i<LABELS.length; i++){
					if(parser.termEquals(LABELS[i], 0)){
						require(!seen[i], "A member length bin must occur exactly once");
						seen[i]=true; assigned[i]=count;
					}
				}
			}
		}finally{require(!members.close(), "Member census read failed");}
		long all=0, retained=0;
		final ByteBuilder table=new ByteBuilder("length_bin\tpopulation\tassigned\tresidual\tassigned_fraction\n");
		for(int i=0; i<total.length; i++){
			require(seen[i] && total[i]>0 && assigned[i]<=total[i], "Assigned members must fit each populated length bin");
			all=Math.addExact(all, total[i]); retained=Math.addExact(retained, assigned[i]);
			table.append(LABELS[i].substring(8)).tab().append(total[i]).tab().append(assigned[i]).tab()
				.append(total[i]-assigned[i]).tab().append(assigned[i]/(double)total[i], 9).nl();
		}
		require(all==expected && retained==recordedMembers, "Coverage bins must conserve both population and member totals");
		Files.createDirectory(out);
		final ByteStreamWriter writer=new ByteStreamWriter(out.resolve("coverage.tsv").toString(), false, false, false);
		writer.start(); writer.print(table);
		require(!writer.poisonAndWait(), "Coverage table write failed");
		Files.write(out.resolve("PASS"), new ByteBuilder("FAMILY_LENGTH_COVERAGE_PASS population=").append(all)
			.append(" assigned=").append(retained).append(" residual=").append(all-retained).nl().toBytes());
		System.out.print(table.toString());
	}

	/** Reads each canonical bucket once and verifies totals before publishing small tables. */
	private static void count(final Path config, final String pin, final Path out, final VerificationMode verification) throws Exception{
		final BasicFileAttributes initial;
		if(verification==VerificationMode.STRUCTURAL){
			require(pin.equals("unmeasured") || pin.matches("[0-9a-f]{20}"), "Invalid structural preparation identity: "+pin);
			initial=structuralState(config);
		}else{
			require(pin.matches("[0-9a-f]{20}") && pin.equals(DigestSuffix.file(config.toString())), "Preparation config fingerprint differs");
			initial=null;
		}
		require(!Files.exists(out), "Census output must be a fresh directory");
		final HashMap<String,String> settings=readConfig(config);
		require(required(settings, "schema").equals("family_preparation_v1"), "Expected the native family preparation schema");
		Shared.AMINO_IN=true;
		Shared.setThreads(2);
		ReadWrite.USE_UNPIGZ=false;
		final int expectedFamilies=Integer.parseInt(required(settings, "expectedfamilies"));
		final long expectedMembers=Long.parseLong(required(settings, "expectedmembers"));
		final MemberLengthRange range=MemberLengthRange.from(settings);
		final String[] labels=range.min<100 ? EXTENDED_LABELS : LABELS;
		final MemberInputs source=new MemberInputs(memberPath(config, settings, "families"), required(settings, "familiessha80"),
			memberPath(config, settings, "chunks"), required(settings, "chunkssha80"),
			memberPath(config, settings, "consensus"), required(settings, "consensussha80"), expectedFamilies, expectedMembers,
			Integer.parseInt(required(settings, "buckets")), Integer.parseInt(required(settings, "batches")),
			Integer.parseInt(required(settings, "maxreferenceaa")), range, verification);
		final long maxMembers=Long.parseLong(required(settings, "maxbucketmembers"));
		final long maxBytes=Long.parseLong(required(settings, "maxbucketbytes"));
		final long[] bins=new long[labels.length];
		final int[] lengths=new int[201];
		final ByteBuilder familyTable=new ByteBuilder("family\tmembers\tlength_p10\tlength_median\tlength_p90\tcount_under_200\n");
		final ByteBuilder meanTable=new ByteBuilder("family\tmembers\tlength_sum\tlength_mean\n");
		final boolean[] seen=new boolean[expectedFamilies];
		long families=0, members=0, fiveShort=0;
		for(int bucket=0; bucket<source.buckets; bucket++){
			final MemberInput[] rows=source.loadBucket(bucket, maxMembers, maxBytes);
			for(int i=0; i<rows.length; i++){
				final MemberInput row=rows[i];
				require(row.rank>=0 && row.rank<seen.length && !seen[row.rank], "Canonical physical family ranks must be counted once");
				seen[row.rank]=true;
				java.util.Arrays.fill(lengths, 0);
				final int shortCount=addFamily(row.members, bins, lengths, range);
				final long lengthSum=lengthSum(lengths, row.members.length, range);
				meanTable.append(row.id).tab().append(row.members.length).tab().append(lengthSum).tab();
				if(row.members.length==0){meanTable.append("NA");}
				else{meanTable.append(lengthSum/(double)row.members.length, 9);}
				meanTable.nl();
				if(shortCount>=5){fiveShort++;}
				familyTable.append(row.id).tab().append(row.members.length).tab();
				if(row.members.length==0){familyTable.append("NA\tNA\tNA");}
				else{
					familyTable.append(quantile(lengths, row.members.length, 10, range), 1).tab()
						.append(quantile(lengths, row.members.length, 50, range), 1).tab()
						.append(quantile(lengths, row.members.length, 90, range), 1);
				}
				familyTable.tab().append(shortCount).nl();
				members=Math.addExact(members, row.members.length);
				families++;
				rows[i]=null;
			}
		}
		source.verifyComplete();
		require(families==expectedFamilies && members==expectedMembers, "Member census must conserve the pinned roster totals");
		long binned=0;
		for(long value:bins){binned=Math.addExact(binned, value);}
		require(binned==members && fiveShort<=families, "Length bins and qualifying family count must partition the input");
		if(initial==null){require(pin.equals(DigestSuffix.file(config.toString())), "Preparation config changed during census");}
		else{verifyStable(config, initial);}
		Files.createDirectory(out);
		final ByteBuilder table=new ByteBuilder("metric\tcount\n");
		for(int i=0; i<bins.length; i++){table.append(labels[i]).tab().append(bins[i]).nl();}
		table.append("families_with_at_least_5_members_below_200").tab().append(fiveShort).nl();
		table.append("families").tab().append(families).nl().append("members").tab().append(members).nl();
		final ByteStreamWriter writer=new ByteStreamWriter(out.resolve("summary.tsv").toString(), false, false, false);
		writer.start(); writer.print(table);
		require(!writer.poisonAndWait(), "Member census output failed");
		final ByteStreamWriter familyWriter=new ByteStreamWriter(out.resolve("family_lengths.tsv").toString(), false, false, false);
		familyWriter.start(); familyWriter.print(familyTable);
		require(!familyWriter.poisonAndWait(), "Per-family member lengths write failed");
		final ByteStreamWriter meanWriter=new ByteStreamWriter(out.resolve("family_mean_lengths.tsv").toString(), false, false, false);
		meanWriter.start(); meanWriter.print(meanTable);
		require(!meanWriter.poisonAndWait(), "Per-family member length sums/means write failed");
		Files.write(out.resolve("quantiles.tsv"), ("p10_p90\tnearest_rank\nmedian\tmean_of_two_middle_members_if_even\n"
			+"empty_family\tNA_quantiles_with_zero_members\nlength_source\tretained_member_peptides_not_consensus_or_original_gene\n"+range.recipe("\t"))
			.getBytes(StandardCharsets.US_ASCII));
		Files.write(out.resolve("PASS"), new ByteBuilder("FAMILY_MEMBER_LENGTH_CENSUS_PASS families=").append(families)
			.append(" members=").append(members).nl().toBytes());
		System.out.print(table.toString());
	}

	/** Adds one whole family to primitive length totals and returns its count below200aa. */
	public static int addFamily(final byte[][] members, final long[] bins){
		return addFamily(members, bins, null);
	}

	/** Optionally fills a caller-owned histogram in the same traversal as the global census. */
	public static int addFamily(final byte[][] members, final long[] bins, final int[] histogram){
		return addFamily(members, bins, histogram, MemberLengthRange.LEGACY);
	}

	/** Counts declared-range members; ranges below100 add three bins ahead of the legacy four. */
	public static int addFamily(final byte[][] members, final long[] bins, final int[] histogram, final MemberLengthRange range){
		require(range!=null && members!=null && bins!=null && bins.length==(range.min<100 ? 7 : 4),
			"Family counters must match the declared member range and its disjoint length bins");
		require(histogram==null || histogram.length==201, "Member histogram must address every allowed40..200 length");
		int shortCount=0;
		for(byte[] member:members){
			require(member!=null, "Canonical MemberInputs must provide every declared member");
			final int bin=lengthBin(member.length, range);
			bins[bin]=Math.addExact(bins[bin], 1);
			if(histogram!=null){histogram[member.length]=Math.addExact(histogram[member.length], 1);}
			if(member.length<200){shortCount++;}
		}
		assert(shortCount<=members.length) : "Each retained member contributes at most once to a family's short count";
		return shortCount;
	}

	/** Exact retained-residue sum, paired with count for unrounded family-mean decisions.
	 * The printed mean is descriptive; consumers compare sum with boundary*count.
	 * Empty families have sum/count0/0 and an explicit unavailable mean, never an invented zero. */
	public static long lengthSum(final int[] histogram, final int count, final MemberLengthRange range){
		require(range!=null && histogram!=null && histogram.length==201 && count>=0,
			"Member length sums need a complete declared-range histogram and nonnegative count");
		long observed=0, residues=0;
		for(int length=0; length<histogram.length; length++){
			final int n=histogram[length];
			require(n>=0 && (range.contains(length) || n==0),
				"Family mean cannot include negative counts or lengths outside its declared retained-member range");
			observed=Math.addExact(observed, n);
			residues=Math.addExact(residues, (long)length*n);
		}
		require(observed==count, "Family mean denominator must equal the complete retained-member population");
		assert(count>0 || residues==0) : "An empty validated histogram cannot contribute residues to a family mean";
		return residues;
	}

	/** Nearest-rank deciles; the median averages the two central observations for even populations. */
	public static double quantile(final int[] histogram, final int count, final int percent){
		return quantile(histogram, count, percent, MemberLengthRange.LEGACY);
	}

	/** Computes quantiles only from lengths inside the declared range, preserving the legacy convention. */
	public static double quantile(final int[] histogram, final int count, final int percent, final MemberLengthRange range){
		require(range!=null && histogram!=null && histogram.length==201 && count>0 && (percent==10 || percent==50 || percent==90),
			"Quantiles require a populated declared-range histogram and a requested10/50/90 percentile");
		long sum=0;
		for(int i=0; i<histogram.length; i++){
			require(histogram[i]>=0 && (range.contains(i) || histogram[i]==0), "Length histogram contains a count outside the declared member range");
			sum+=histogram[i];
		}
		require(sum==count, "Quantile denominator must equal the complete family member population");
		if(percent==50){return (atRank(histogram, (count+1L)/2, range)+atRank(histogram, (count+2L)/2, range))/2.0;}
		return atRank(histogram, ((long)count*percent+99)/100, range);
	}

	/** Returns the length at a one-based order statistic without sorting member sequences. */
	private static int atRank(final int[] histogram, final long rank, final MemberLengthRange range){
		require(rank>0, "Member order statistics are one-based");
		long cumulative=0;
		for(int length=range.min; length<=range.max; length++){
			cumulative+=histogram[length];
			if(cumulative>=rank){return length;}
		}
		throw new IllegalArgumentException("Member order statistic exceeds histogram population: "+rank);
	}

	/** The native member reader contract is100..200aa; consensus length is never used here. */
	public static int lengthBin(final int length){
		return lengthBin(length, MemberLengthRange.LEGACY);
	}

	/** Explicit short ranges use40..59,60..79,80..99 followed by the unchanged legacy bins. */
	public static int lengthBin(final int length, final MemberLengthRange range){
		require(range!=null && range.contains(length), "Member length must fit its declared reader range: "+length+" in "+range);
		if(length<100){return length<60 ? 0 : length<80 ? 1 : 2;}
		return (range.min<100 ? 3 : 0)+(length<125 ? 0 : length<150 ? 1 : length<200 ? 2 : 3);
	}

	/** Reads the small pinned native config; unrelated preparation options remain untouched. */
	private static HashMap<String,String> readConfig(final Path path){
		final HashMap<String,String> settings=new HashMap<String,String>();
		final ByteFile input=ByteFile.makeByteFile1(path.toString(), false);
		try{
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				final String text=new String(line, StandardCharsets.UTF_8).trim();
				if(text.isEmpty() || text.startsWith("#")){continue;}
				final String[] pair=Parse.splitOnFirst(text, '=');
				require(pair.length==2 && !pair[1].isEmpty() && settings.put(pair[0], pair[1])==null, "Invalid or duplicate native config key");
			}
		}finally{require(!input.close(), "Native config read failed");}
		return settings;
	}

	/** Missing metadata must never be replaced with guessed counts or paths. */
	private static String required(final HashMap<String,String> values, final String key){
		final String value=values.get(key);
		require(value!=null && !value.isEmpty(), "Missing required field: "+key);
		return value;
	}
	/** Native preparation resolves relative member manifests against the configuration directory. */
	private static Path memberPath(final Path config, final HashMap<String,String> values, final String key){
		final Path path=Paths.get(required(values, key));
		return path.isAbsolute() ? path : config.toAbsolutePath().getParent().resolve(path).normalize();
	}
	/** Structural checks are the default; hashing requires an explicit exact-match request. */
	private static VerificationMode verification(final HashMap<String,String> options, final String mode){
		final String policy=options.getOrDefault("verification", "structural");
		require(policy.equals("sha80") || policy.equals("structural"),
			"Unknown "+mode+" census verification policy: "+policy);
		return policy.equals("sha80") ? VerificationMode.SHA80 : VerificationMode.STRUCTURAL;
	}
	/** Verifies metadata after reading; this does not establish byte identity. */
	private static void verifyStable(final Path path, final BasicFileAttributes initial) throws java.io.IOException{
		final BasicFileAttributes current=structuralState(path);
		require(initial.size()==current.size() && initial.lastModifiedTime().equals(current.lastModifiedTime()) &&
			java.util.Objects.equals(initial.fileKey(), current.fileKey()), "Census input metadata changed while reading: "+path);
	}
	/** Structural mode binds the regular input's identity and size, without claiming byte identity. */
	private static BasicFileAttributes structuralState(final Path path) throws java.io.IOException{
		require(path.isAbsolute() && !Files.isSymbolicLink(path), "Structural census requires an absolute nonsymlink input: "+path);
		final BasicFileAttributes state=Files.readAttributes(path, BasicFileAttributes.class, LinkOption.NOFOLLOW_LINKS);
		require(state.isRegularFile(), "Structural census input must be a regular file: "+path);
		return state;
	}
	private static void require(final boolean valid, final String reason){if(!valid){throw new IllegalArgumentException(reason);}}
	private static final String[] LABELS={"members_100_124", "members_125_149", "members_150_199", "members_200"};
	private static final String[] EXTENDED_LABELS={"members_40_59", "members_60_79", "members_80_99",
		"members_100_124", "members_125_149", "members_150_199", "members_200"};
}
