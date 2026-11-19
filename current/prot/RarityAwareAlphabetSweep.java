package prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.LinkOption;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.attribute.BasicFileAttributes;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Locale;
import java.util.Objects;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Sweeps complete-link alphabets with a frequency-uniformity merge term.
 *
 * @author Ayaka
 */
public final class RarityAwareAlphabetSweep{

	private RarityAwareAlphabetSweep(){}

	public static void main(String[] args) throws Exception{
		final PreParser pp=new PreParser(args, RarityAwareAlphabetSweep.class, false);
		try{run(pp.args);}
		finally{Shared.closeStream(pp.outstream);}
	}

	public static void run(String[] args) throws Exception{
		final HashMap<String,String> options=parse(args);
		final Path input=Paths.get(required(options, "in")).toAbsolutePath().normalize();
		final Path output=Paths.get(required(options, "out")).toAbsolutePath().normalize();
		final String pin=options.get("insha80");
		final String diagnosticFlag=options.getOrDefault("mergediagnostics", "false");
		require(diagnosticFlag.equals("true") || diagnosticFlag.equals("false"), "mergediagnostics must be true or false");
		final boolean diagnostics=Boolean.parseBoolean(diagnosticFlag);
		final String verification=options.containsKey("verification") ? options.get("verification") :
			"normalized_probabilities".equals(options.get("format")) ? "structural" : "sha80";
		require(verification.equals("sha80") || verification.equals("structural"),
			"verification must be sha80 or structural");
		final boolean structural=verification.equals("structural");
		require(!structural || pin==null, "Structural verification does not accept insha80; retain inherited identity externally");
		if(pin!=null){
			DigestSuffix.requireSuffix(pin, "Dimer substitution table");
			require(pin.equals(DigestSuffix.file(input.toString())), "Dimer substitution table differs from its declared sha80");
		}
		require(Files.isRegularFile(input) && !Files.isSymbolicLink(input), "Input must be a regular dimer table: "+input);
		require(!Files.exists(output) && Files.isDirectory(output.getParent()),
			"Output directory must be fresh and under an existing parent: "+output);
		final InputState before=structural ? new InputState(input) : null;
		final double[] lambdas=doubles(options.containsKey("lambdas") ? options.get("lambdas") : "0,0.1,0.25,0.5,1");
		final int[] targets=ints(options.containsKey("targets") ? options.get("targets") : "120");
		final String format=options.containsKey("format") ? options.get("format") : "symmetric";
		require(format.equals("symmetric") || format.equals("member_pair_counts") || format.equals("member_pair_blend")
			|| format.equals("member_pair_band_blend") || format.equals("normalized_probabilities"),
			"format must be symmetric, member_pair_counts, member_pair_blend, member_pair_band_blend or normalized_probabilities");
		final String band=options.containsKey("band") ? options.get("band") : "<40";
		final Table table=format.equals("symmetric") ? read(input,band) :
			format.equals("member_pair_counts") ? readMemberPairs(input,band) :
			format.equals("member_pair_blend") ? readMemberBlend(input,band) :
			format.equals("normalized_probabilities") ? readNormalizedProbabilities(input,band) : readMemberBandBlend(input,band);
		for(double lambda : lambdas){require(lambda>=0 && Double.isFinite(lambda), "Lambda must be finite and nonnegative: "+lambda);}
		for(int target : targets){require(target>0 && target<=table.n, "Target must be in [1,"+table.n+"]: "+target);}
		Files.createDirectory(output);
		write(output.resolve("alphabets.tsv"), table, lambdas, targets, diagnostics);
		final ByteBuilder manifest=new ByteBuilder("key\tvalue\nsource\t").append(input.toString()).nl();
		if(pin!=null){manifest.append("source_sha80\t").append(pin).nl();}
		manifest.append("schema\trarity_aware_alphabet_sweep_v1\n")
			.append("input_format\t").append(format).nl()
			.append("band\t").append(table.band).nl()
			.append("symbols\t").append(table.n).nl()
			.append("substitution_similarity\tcomplete-link weakest original odds / best_candidate_weakest_original_odds\n")
			.append("uniformity\t-sum((group_frequency - 1/group_count)^2)\n")
			.append("frequency_uniformity_gain\tpost_merge_uniformity - pre_merge_uniformity\n")
			.append("normalized_frequency_uniformity_gain\tfrequency_uniformity_gain / max_abs_candidate_gain_per_step\n")
			.append("merge_score\tnormalized_substitution_similarity + lambda * normalized_frequency_uniformity_gain\n")
			.append("tie_breaker\texisting FamilyDimerSubstitutionTable whole-partition lexical order\n")
			.append("lambda0_contract\texact FamilyDimerSubstitutionTable.reduce replay\n");
		if(format.equals("member_pair_counts")){manifest.append("input_view\tpooled\n");}
		if(format.equals("member_pair_blend") || format.equals("member_pair_band_blend")){manifest.append("input_view\tpooled\n")
			.append("endpoint_frequency\texact_double_marginals_from_equal_source_probability_blend;no_integer_rounding\n");}
		if(format.equals("normalized_probabilities")){manifest.append("input_view\tpooled\n")
			.append("endpoint_frequency\tdouble_marginals_from_recursive_normalized_counts;no_integer_rounding\n");}
		if(structural){
			manifest.append("verification\tstructural\n")
				.append("input_bytes\t").append(before.attributes.size()).nl()
				.append("input_file_key\t").append(String.valueOf(before.attributes.fileKey())).nl()
				.append("input_modified\t").append(before.attributes.lastModifiedTime().toString()).nl()
				.append("input_created\t").append(before.attributes.creationTime().toString()).nl();
		}else{
			manifest.append("alphabets_sha80\t").append(DigestSuffix.file(output.resolve("alphabets.tsv").toString())).nl();
		}
		write(output.resolve("manifest.tsv"), manifest);
		if(pin!=null){require(pin.equals(DigestSuffix.file(input.toString())), "Dimer substitution table changed during sweep");}
		if(structural){require(before.matches(new InputState(input)), "Dimer substitution table metadata changed during sweep");}
		write(output.resolve("PASS"), new ByteBuilder("RARITY_AWARE_ALPHABET_SWEEP_PASS symbols=")
			.append(table.n).append(" lambdas=").append(lambdas.length).append(" targets=").append(targets.length).nl());
	}

	/** Loads explicitly versioned recursive probability state, preserving fractional endpoint marginals. */
	private static Table readNormalizedProbabilities(Path path, String band){
		final FamilySubstitutionCountBlend.Table source=FamilySubstitutionCountBlend.readProbabilities(path, true, new String[]{band}).get(band);
		final Table table=new Table(band);
		for(String label:FamilyDimerSubstitutionTable.dimerLabels()){table.id(label);}
		table.n=source.n;
		table.separator=",";
		table.odds=source.odds();
		table.frequency=source.endpointFrequencies();
		table.verifyLabels(table.labels());
		return table;
	}

	/** Read the complete pooled ordered matrix without recounting or an intermediate table. */
	private static Table readMemberPairs(Path path, String band){
		require(band!=null && band.length()>0, "Band is required");
		final String[] labels=FamilyDimerSubstitutionTable.dimerLabels();
		final int n=labels.length;
		final long[] ordered=new long[n*n];
		final LineParser1 parser=new LineParser1('\t');
		final ByteFile reader=ByteFile.makeByteFile(path.toString(), false);
		int cells=0;
		try{
			final byte[] header=reader.nextLine();
			require(header!=null && new String(header, StandardCharsets.US_ASCII).equals(PAIR_HEADER),
				"Expected FamilyMemberPairDimerExport ordered-count schema");
			for(byte[] line=reader.nextLine(); line!=null; line=reader.nextLine()){
				parser.set(line);
				require(parser.terms()==5, "Member-pair count row requires five fields");
				if(!parser.termEquals("pooled", 0) || !parser.termEquals(band, 1)){continue;}
				require(cells<ordered.length && parser.termEquals(labels[cells/n], 2)
					&& parser.termEquals(labels[cells%n], 3),
					"Pooled ordered counts must cover each canonical row-major cell once; cell="+cells);
				final long count=parser.parseLong(4);
				require(count>=0, "Dimer observations cannot be negative; cell="+cells);
				ordered[cells++]=count;
			}
		}finally{require(!reader.close(), "Failed reading member-pair dimer counts");}
		require(cells==ordered.length, "Incomplete pooled dimer matrix in band "+band+": "+cells);
		// Reuse the original endpoint symmetrization and double arithmetic verbatim.
		final FamilyDimerSubstitutionTable.Table original=new FamilyDimerSubstitutionTable.Table(ordered, n);
		final Table table=new Table(band);
		for(String label : labels){table.id(label);}
		for(int i=0; i<n; i++){
			for(int j=i; j<n; j++){
				table.rows.add(new Cell(i, j, original.endpoints[i], original.endpoints[j], original.odds[i*n+j]));
			}
		}
		table.finish();
		return table;
	}

	/** Read the accepted equal-source probability blend without inventing integer endpoint counts. */
	private static Table readMemberBlend(Path path,String band){
		require("<40".equals(band),"The retained blend-frequency schema binds the original below40 population");
		final String[] labels=FamilyDimerSubstitutionTable.dimerLabels();final int n=labels.length;
		final double[] mixed=new double[n*n],savedOdds=new double[n*n],frequency=new double[n];
		final ByteFile reader=ByteFile.makeByteFile(path.toString(),false);final LineParser1 parser=new LineParser1('\t');
		double firstTotal=0,secondTotal=0;
		try{
			require(Arrays.equals(reader.nextLine(),BLEND_HEADER.getBytes(StandardCharsets.US_ASCII)),
				"Require the original six-column consensus/pairwise/blend frequency table");
			for(int i=0;i<n;i++){for(int j=i;j<n;j++){
				final byte[] line=reader.nextLine();require(line!=null,"The blend must contain every canonical upper-triangle cell");parser.set(line);
				require(parser.terms()==6 && parser.termEquals(labels[i],0) && parser.termEquals(labels[j],1),
					"Blend cells must occur exactly once in canonical upper-triangle order");
				final double first=Double.parseDouble(parser.parseString(2)),second=Double.parseDouble(parser.parseString(3)),
					value=Double.parseDouble(parser.parseString(4)),odds=Double.parseDouble(parser.parseString(5));
				require(Double.isFinite(first) && Double.isFinite(second) && Double.isFinite(value) && Double.isFinite(odds)
					&& first>=0 && second>=0 && value>=0 && odds>=0,"Blend probabilities and odds must remain finite and nonnegative");
				require(Double.doubleToLongBits(value)==Double.doubleToLongBits(.5*first+.5*second),
					"The released blend is exactly half of each separately normalized source, not a count-weighted mixture");
				firstTotal+=first;secondTotal+=second;mixed[i*n+j]=value;savedOdds[i*n+j]=savedOdds[j*n+i]=odds;
				final double half=.5*value;frequency[i]+=half;frequency[j]+=half;
			}}
			require(reader.nextLine()==null,"Unexpected trailing blend-frequency row");
		}finally{require(!reader.close(),"Retained blend-frequency input failed");}
		require(Math.abs(firstTotal-1)<1e-10 && Math.abs(secondTotal-1)<1e-10,
			"Both retained sources must have unit probability mass before blending");
		final double[] expected=FamilyMemberPairDimerExport.odds(mixed);
		for(int i=0;i<expected.length;i++){
			require(Double.doubleToLongBits(expected[i])==Double.doubleToLongBits(savedOdds[i]),
				"Blend odds must exactly reproduce the accepted endpoint arithmetic; cell="+i);
		}
		// These are probabilities already; finish() deliberately retains its integer contract for legacy readers.
		final Table table=new Table(band);for(String label:labels){table.id(label);}
		table.n=n;table.separator=",";table.odds=savedOdds;table.frequency=frequency;table.verifyLabels(table.labels());
		double total=0;for(double value:frequency){require(Double.isFinite(value) && value>=0,"Every mixed endpoint must be a valid probability");total+=value;}
		require(Math.abs(total-1)<1e-10,"Mixed endpoint probabilities must sum to one without integer rescaling");
		return table;
	}

	/** Versioned band-specific blend; retain the original six-column reader and its population guard. */
	private static Table readMemberBandBlend(Path path, String band){
		require(band!=null && band.length()>0, "Band is required for a band-specific blend");
		final String[] labels=FamilyDimerSubstitutionTable.dimerLabels();
		final int n=labels.length;
		final double[] mixed=new double[n*n];
		final double[] savedOdds=new double[n*n];
		final double[] frequency=new double[n];
		final ByteFile reader=ByteFile.makeByteFile(path.toString(), false);
		final LineParser1 parser=new LineParser1('\t');
		double firstTotal=0;
		double secondTotal=0;
		try{
			require(Arrays.equals(reader.nextLine(), BAND_BLEND_MAGIC.getBytes(StandardCharsets.US_ASCII)),
				"Require bbgap_dimer_blend version1 before interpreting band-specific probabilities");
			require(Arrays.equals(reader.nextLine(), BAND_BLEND_HEADER.getBytes(StandardCharsets.US_ASCII)),
				"Require the seven-column band/consensus/pairwise/blend frequency table");
			for(int i=0; i<n; i++){
				for(int j=i; j<n; j++){
					final byte[] line=reader.nextLine();
					require(line!=null, "The band blend must contain every canonical upper-triangle cell");
					parser.set(line);
					require(parser.terms()==7 && parser.termEquals(band, 0)
						&& parser.termEquals(labels[i], 1) && parser.termEquals(labels[j], 2),
						"Band blend cells must use the requested band and canonical upper-triangle order; band="+band+", cell="+i+","+j);
					final double first=Double.parseDouble(parser.parseString(3));
					final double second=Double.parseDouble(parser.parseString(4));
					final double value=Double.parseDouble(parser.parseString(5));
					final double odds=Double.parseDouble(parser.parseString(6));
					require(Double.isFinite(first) && Double.isFinite(second) && Double.isFinite(value) && Double.isFinite(odds)
						&& first>=0 && second>=0 && value>=0 && odds>=0,
						"Blend probabilities and odds must remain finite and nonnegative");
					require(Double.doubleToLongBits(value)==Double.doubleToLongBits(.5*first+.5*second),
						"The band blend is exactly half of each separately normalized source, not a count-weighted mixture");
					firstTotal+=first;
					secondTotal+=second;
					mixed[i*n+j]=value;
					savedOdds[i*n+j]=odds;
					savedOdds[j*n+i]=odds;
					final double half=.5*value;
					frequency[i]+=half;
					frequency[j]+=half;
				}
			}
			require(reader.nextLine()==null, "Unexpected trailing band-blend-frequency row");
		}finally{
			require(!reader.close(), "Band-blend-frequency input failed");
		}
		require(Math.abs(firstTotal-1)<1e-10 && Math.abs(secondTotal-1)<1e-10,
			"Both band-specific sources must have unit probability mass before blending");
		final double[] expected=FamilyMemberPairDimerExport.odds(mixed);
		for(int i=0; i<expected.length; i++){
			require(Double.doubleToLongBits(expected[i])==Double.doubleToLongBits(savedOdds[i]),
				"Band blend odds must exactly reproduce the accepted endpoint arithmetic; cell="+i);
		}
		final Table table=new Table(band);
		for(String label : labels){table.id(label);}
		table.n=n;
		table.separator=",";
		table.odds=savedOdds;
		table.frequency=frequency;
		table.verifyLabels(table.labels());
		double total=0;
		for(double value : frequency){
			require(Double.isFinite(value) && value>=0, "Every mixed endpoint must be a valid probability");
			total+=value;
		}
		require(Math.abs(total-1)<1e-10, "Mixed endpoint probabilities must sum to one without integer rescaling");
		return table;
	}

	private static Table read(Path path, String band){
		require(band!=null && band.length()>0, "Band is required");
		final Table table=new Table(band);
		final LineParser1 parser=new LineParser1('\t');
		final ByteFile reader=ByteFile.makeByteFile(path.toString(), false);
		try{
			final byte[] header=reader.nextLine();
			require(header!=null && new String(header, StandardCharsets.US_ASCII).equals(HEADER),
				"Expected exact retained dimer-table schema");
			for(byte[] line=reader.nextLine(); line!=null; line=reader.nextLine()){
				parser.set(line); require(parser.terms()==9, "Dimer-table row requires nine fields");
				if(!parser.termEquals(band, 0)){continue;}
				final String left=parser.parseString(1), right=parser.parseString(2);
				final int i=table.id(left), j=table.id(right);
				require(i<=j, "Dimer table must contain upper-triangle cells: "+left+","+right);
				// The accepted retained table writes Double.toString odds; parse them exactly to preserve saved ties.
				table.rows.add(new Cell(i, j, parser.parseLong(4), parser.parseLong(5),
					Double.parseDouble(parser.parseString(8))));
			}
		}finally{require(!reader.close(), "Failed reading dimer substitution table");}
		table.finish();
		return table;
	}

	private static void write(Path path, Table table, double[] lambdas, int[] targets, boolean diagnostics){
		final Diagnostics trace=diagnostics ? new Diagnostics(table) : null;
		final ByteStreamWriter writer=new ByteStreamWriter(path.toString(), false, false, false);
		final ByteBuilder line=new ByteBuilder();
		writer.start();
		try{
			writer.println("band\tlambda\ttarget\tgroups");
			final FamilyDimerSubstitutionTable.Hierarchy lambda0=contains(lambdas, 0) ?
				FamilyDimerSubstitutionTable.reduce(table.odds, table.labels(), table.separator(), targets) : null;
			for(double lambda : lambdas){
				final Hierarchy current=reduce(table, lambda, trace, targets);
				if(lambda==0){require(same(lambda0, current), "Lambda 0 diverged from FamilyDimerSubstitutionTable.reduce");}
				for(Cut cut : current.cuts){
					line.clear().append(table.band).tab().append(Double.toString(lambda)).tab().append(cut.symbols).tab()
						.append(join(cut.groups));
					writer.println(line);
				}
			}
		}finally{require(!writer.poisonAndWait(), "Rarity-aware alphabet output failed: "+path);}
		if(trace!=null){trace.write(path.getParent());}
	}

	private static Hierarchy reduce(Table table, double lambda, Diagnostics trace, int... targets){
		final int n=table.n;
		final Group[] groups=new Group[n];
		final String[] labels=table.labels();
		table.verifyLabels(labels);
		final String[] mergedCache=new String[n*n];
		final double[] scores=table.odds.clone();
		int minimum=n;
		for(int target : targets){minimum=Math.min(minimum, target);}
		for(int i=0; i<n; i++){groups[i]=new Group(new int[]{i}, labels[i], table.frequency[i]);}
		final Hierarchy result=new Hierarchy();
		final PartitionCursor left=new PartitionCursor(), right=new PartitionCursor();
		int active=n;
		if(trace!=null){trace.cut(groups, targets, active, lambda);}
		capture(result, groups, targets, active);
		while(active>minimum){
			final String[] sorted=texts(groups);
			final boolean raw=lambda==0;
			final double sum=sum(groups), sumSquares=sumSquares(groups), base=uniformity(sum, sumSquares, active),
				oddsScale=raw ? 1 : oddsScale(scores, groups), gainScale=raw ? 0 : gainScale(groups, active, sum, sumSquares, base);
			double best=Double.NEGATIVE_INFINITY;
			int bestA=-1, bestB=-1, maxTies=0;
			for(int a=0; a<n; a++){
				if(groups[a]==null){continue;}
				for(int b=a+1; b<n; b++){
					if(groups[b]==null){continue;}
					final double value;
					if(raw){value=scores[a*n+b];}
					else{
						final double rawGain=uniformityAfter(groups, active, sum, sumSquares, a, b)-base;
						final double gain=gainScale>0 ? rawGain/gainScale : 0;
						value=scores[a*n+b]/oddsScale+lambda*gain;
					}
					if(value<best){continue;}
					if(trace!=null){maxTies=value>best ? 1 : maxTies+1;}
					if(value>best || bestA<0 || before(groups, labels, table.separator, mergedCache, sorted, a, b, bestA, bestB, left, right)){
						best=value; bestA=a; bestB=b;
					}
				}
			}
			assert(bestA>=0 && bestB>bestA) : "At least two active groups supply a merge candidate";
			if(trace!=null){trace.merge(groups, scores, lambda, n-active+1, active, bestA, bestB, best,
				sum, sumSquares, base, oddsScale, gainScale, maxTies);}
			final Group a=groups[bestA], b=groups[bestB];
			for(int c=0; c<n; c++){
				if(c!=bestA && c!=bestB && groups[c]!=null){
					scores[bestA*n+c]=scores[c*n+bestA]=Math.min(scores[bestA*n+c], scores[bestB*n+c]);
				}
				mergedCache[bestA*n+c]=mergedCache[c*n+bestA]=null;
				mergedCache[bestB*n+c]=mergedCache[c*n+bestB]=null;
			}
			groups[bestA]=new Group(union(a.members, b.members), mergedText(groups, labels, table.separator, mergedCache, bestA, bestB), a.frequency+b.frequency);
			groups[bestB]=null;
			active--;
			if(trace!=null){trace.cut(groups, targets, active, lambda);}
			capture(result, groups, targets, active);
		}
		return result;
	}

	private static double oddsScale(double[] scores, Group[] groups){
		final int n=groups.length;
		double best=Double.NEGATIVE_INFINITY;
		for(int a=0; a<n; a++){
			if(groups[a]==null){continue;}
			for(int b=a+1; b<n; b++){
				if(groups[b]!=null){best=Math.max(best, scores[a*n+b]);}
			}
		}
		require(best>=0 && Double.isFinite(best), "Every merge step needs a finite nonnegative best odds");
		return best>0 ? best : 1;
	}

	private static double gainScale(Group[] groups, int active, double sum, double sumSquares, double base){
		double max=0;
		for(int a=0; a<groups.length; a++){
			if(groups[a]==null){continue;}
			for(int b=a+1; b<groups.length; b++){
				if(groups[b]!=null){max=Math.max(max, Math.abs(uniformityAfter(groups, active, sum, sumSquares, a, b)-base));}
			}
		}
		return max;
	}

	private static double uniformityAfter(Group[] groups, int active, double sum, double sumSquares, int a, int b){
		final double nextTarget=1.0/(active-1);
		final double af=groups[a].frequency, bf=groups[b].frequency, mf=af+bf;
		return -((sumSquares-af*af-bf*bf+mf*mf)-2*nextTarget*sum+(active-1)*nextTarget*nextTarget);
	}

	private static double uniformity(double sum, double sumSquares, int active){
		final double target=1.0/active;
		return -(sumSquares-2*target*sum+active*target*target);
	}

	private static double sum(Group[] groups){
		double sum=0;
		for(Group group : groups){
			if(group!=null){sum+=group.frequency;}
		}
		return sum;
	}

	private static double sumSquares(Group[] groups){
		double sum=0;
		for(Group group : groups){
			if(group!=null){sum+=group.frequency*group.frequency;}
		}
		return sum;
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

	private static void capture(Hierarchy result, Group[] groups, int[] targets, int active){
		boolean requested=false;
		for(int target : targets){if(target==active){requested=true;}}
		if(!requested){return;}
		result.cuts.add(new Cut(active, texts(groups)));
	}

	private static boolean same(FamilyDimerSubstitutionTable.Hierarchy old, Hierarchy current){
		if(old.cuts.size()!=current.cuts.size()){return false;}
		for(int i=0; i<old.cuts.size(); i++){
			if(old.cuts.get(i).symbols!=current.cuts.get(i).symbols || !Arrays.equals(old.cuts.get(i).groups, current.cuts.get(i).groups)){
				return false;
			}
		}
		return true;
	}

	private static boolean contains(double[] array, double value){
		for(double element : array){if(element==value){return true;}}
		return false;
	}

	private static String join(String[] groups){
		final StringBuilder sb=new StringBuilder();
		for(String group : groups){
			if(sb.length()>0){sb.append('/');}
			sb.append(group);
		}
		return sb.toString();
	}

	private static HashMap<String,String> parse(String[] args){
		final HashMap<String,String> map=new HashMap<String,String>();
		for(String arg : args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			require(split.length==2 && split[1].length()>0, "Expected flag=value: "+arg);
			final String key=split[0].toLowerCase(Locale.ROOT);
			require(OPTIONS.contains(key) && map.put(key, split[1])==null, "Unknown or duplicate option: "+key);
		}
		require(map.containsKey("in") && map.containsKey("out"), "Require in and out");
		return map;
	}

	private static double[] doubles(String raw){
		final String[] split=raw.split(",", -1);
		final double[] out=new double[split.length];
		require(out.length>0, "Empty double list");
		for(int i=0; i<out.length; i++){out[i]=Double.parseDouble(split[i]);}
		return out;
	}

	private static int[] ints(String raw){
		final String[] split=raw.split(",", -1);
		final int[] out=new int[split.length];
		require(out.length>0, "Empty int list");
		for(int i=0; i<out.length; i++){out[i]=Integer.parseInt(split[i]);}
		return out;
	}

	private static void write(Path path, ByteBuilder data){
		require(!Files.exists(path), "Refusing to overwrite: "+path);
		final ByteStreamWriter writer=new ByteStreamWriter(path.toString(), false, false, false);
		writer.start(); writer.print(data);
		require(!writer.poisonAndWait(), "Write failed: "+path);
	}

	private static String required(HashMap<String,String> options, String key){
		final String value=options.get(key);
		require(value!=null, "Missing required option: "+key);
		return value;
	}

	private static void require(boolean condition, String message){if(!condition){throw new IllegalArgumentException(message);}}

	/** Observational sidecars; never mutate the live decision cache, groups or scores. */
	private static final class Diagnostics{
		/** Resolve named dimers once and allocate a cache isolated from decision state. */
		Diagnostics(Table table_){
			table=table_;
			labels=table.labels();
			cache=new String[table.n*table.n];
			for(int i=0; i<named.length; i++){
				named[i]=Arrays.asList(labels).indexOf(NAMED[i]);
				require(named[i]>=0, "Diagnostic dimer absent: "+NAMED[i]);
			}
		}
		/** Record a chosen premerge decision and named-group alternatives at its exact scales. */
		void merge(Group[] groups, double[] scores, double lambda, int step, int active, int a, int b,
				double actual, double sum, double squares, double base, double oddsScale, double gainScale, int maxTies){
			Arrays.fill(cache, null);
			final boolean raw=lambda==0;
			final double chosen=score(groups, scores, lambda, active, a, b, sum, squares, base, oddsScale, gainScale);
			require(Double.doubleToLongBits(chosen)==Double.doubleToLongBits(actual), "Diagnostic chosen score drift");
			final String joined=mergedText(groups, labels, table.separator, cache, a, b);
			row(merges, table.band, lambda, step, active, active-1, groups[a].text, groups[b].text, joined,
				scores[a*table.n+b], raw ? "NA" : uniformityAfter(groups, active, sum, squares, a, b)-base,
				raw ? "NA" : oddsScale, raw ? "NA" : gainScale, actual,
				raw ? "raw_odds" : "normalized_odds_plus_lambda_gain", maxTies, "whole_postmerge_partition_lexical");
			final String[] sorted=texts(groups);
			for(int i=0; i<named.length; i++){
				final int owner=contains(groups[a].members, named[i]) ? a : contains(groups[b].members, named[i]) ? b : -1;
				if(owner<0){continue;}
				final int partner=owner==a ? b : a;
				int alternative=-1, ties=0;
				double best=Double.NEGATIVE_INFINITY;
				for(int c=0; c<groups.length; c++){
					if(c==owner || c==partner || groups[c]==null){continue;}
					final int ca=Math.min(owner,c), cb=Math.max(owner,c);
					final double value=score(groups, scores, lambda, active, ca, cb, sum, squares, base, oddsScale, gainScale);
					if(value<best){continue;}
					ties=value>best ? 1 : ties+1;
					if(value>best || alternative<0 || before(groups, labels, table.separator, cache, sorted, ca, cb,
							Math.min(owner,alternative), Math.max(owner,alternative), left, right)){
						best=value;
						alternative=c;
					}
				}
				row(namedMerges, table.band, lambda, step, active, NAMED[i], groups[owner].text, groups[partner].text,
					joined, actual, scores[a*table.n+b], alternative<0 ? "no_alternative" : "available",
					alternative<0 ? "NA" : groups[alternative].text, alternative<0 ? "NA" : best,
					alternative<0 ? "NA" : scores[owner*table.n+alternative], alternative<0 ? "NA" : ties,
					raw ? "raw_odds" : "normalized_odds_plus_lambda_gain", "whole_postmerge_partition_lexical");
			}
		}
		/** Describe singleton-to-group raw odds only at requested cuts, excluding self. */
		void cut(Group[] groups, int[] targets, int active, double lambda){
			boolean requested=false;
			for(int target:targets){
				if(target==active){requested=true;}
			}
			if(!requested){return;}
			for(int i=0; i<named.length; i++){
				Group owner=null;
				for(Group group:groups){
					if(group!=null && contains(group.members,named[i])){
						owner=group;
						break;
					}
				}
				require(owner!=null, "Diagnostic dimer lost from partition");
				final StringBuilder remainder=new StringBuilder();
				double own=Double.POSITIVE_INFINITY;
				for(int member:owner.members){
					if(member!=named[i]){
						if(remainder.length()>0){remainder.append(table.separator);}
						remainder.append(labels[member]);
						own=Math.min(own,table.odds[named[i]*table.n+member]);
					}
				}
				Group alternative=null;
				double best=Double.NEGATIVE_INFINITY;
				for(Group group:groups){
					if(group==null || group==owner){continue;}
					double value=Double.POSITIVE_INFINITY;
					for(int member:group.members){value=Math.min(value,table.odds[named[i]*table.n+member]);}
					if(value>best || (value==best && (alternative==null || group.text.compareTo(alternative.text)<0))){
						best=value;
						alternative=group;
					}
				}
				row(cuts, table.band, lambda, active, NAMED[i], owner.text, owner.members.length==1 ? "singleton" : "available",
					remainder.length()==0 ? "NA" : remainder, owner.members.length==1 ? "NA" : own,
					alternative==null ? "no_other_group" : "available", alternative==null ? "NA" : alternative.text,
					alternative==null ? "NA" : best, "raw_odds_descriptive_not_decision_score", "group_text_lexical");
			}
		}
		/** Replay the actual branch arithmetic without recomputing step normalizers. */
		private double score(Group[] groups, double[] scores, double lambda, int active, int a, int b,
				double sum, double squares, double base, double oddsScale, double gainScale){
			if(lambda==0){return scores[a*table.n+b];}
			final double rawGain=uniformityAfter(groups, active, sum, squares, a, b)-base;
			final double gain=gainScale>0 ? rawGain/gainScale : 0;
			return scores[a*table.n+b]/oddsScale+lambda*gain;
		}
		/** Write only the additive versioned sidecars after alphabet generation succeeds. */
		void write(Path output){
			RarityAwareAlphabetSweep.write(output.resolve("merge_steps.tsv"), new ByteBuilder(merges.toString()));
			RarityAwareAlphabetSweep.write(output.resolve("named_merges.tsv"), new ByteBuilder(namedMerges.toString()));
			RarityAwareAlphabetSweep.write(output.resolve("named_cuts.tsv"), new ByteBuilder(cuts.toString()));
		}
		/** Membership lookup in an immutable, canonically sorted active group. */
		private static boolean contains(int[] members, int id){return Arrays.binarySearch(members,id)>=0;}
		/** Serialize round-trip numeric values and explicit unavailable strings as TSV. */
		private static void row(StringBuilder out, Object... values){
			for(int i=0; i<values.length; i++){
				if(i>0){out.append('\t');}
				out.append(values[i]);
			}
			out.append('\n');
		}
		final Table table;
		final String[] labels, cache;
		final int[] named=new int[NAMED.length];
		final PartitionCursor left=new PartitionCursor(), right=new PartitionCursor();
		final StringBuilder merges=new StringBuilder("#bbgap_alphabet_merge_steps\t1\nband\tlambda\tstep\tactive_before\tactive_after\tleft\tright\tjoined\traw_weakest_odds\traw_frequency_gain\todds_scale\tgain_scale\tactual_score\tscore_units\tmax_score_ties\ttie_breaker\n");
		final StringBuilder namedMerges=new StringBuilder("#bbgap_alphabet_named_merges\t1\nband\tlambda\tstep\tactive_before\tnamed_dimer\tcontaining_group\tchosen_partner\tjoined\tchosen_score\tchosen_raw_odds\talternative_status\talternative_group\talternative_score\talternative_raw_odds\talternative_score_ties\tscore_units\ttie_breaker\n");
		final StringBuilder cuts=new StringBuilder("#bbgap_alphabet_named_cuts\t1\nband\tlambda\ttarget\tnamed_dimer\town_group\town_status\town_other_members\town_raw_odds\talternative_status\talternative_group\talternative_raw_odds\tscore_units\ttie_breaker\n");
		private static final String[] NAMED={"FM","IM","LM","MM","VM","YM","RA","RS","RT"};
	}

	/** Metadata stability brackets for explicit structural verification; no file digest. */
	private static final class InputState{
		/** Capture identity and timestamps, rejecting aliases and non-regular input. */
		InputState(Path path) throws java.io.IOException{
			canonical=path.toRealPath();
			require(canonical.equals(path), "Structural input must have a canonical path: "+path);
			attributes=Files.readAttributes(path, BasicFileAttributes.class, LinkOption.NOFOLLOW_LINKS);
			require(attributes.isRegularFile() && !attributes.isSymbolicLink(), "Structural input must be a regular file: "+path);
		}
		/** Require the same canonical file and metadata at both sweep boundaries. */
		boolean matches(InputState other){
			return canonical.equals(other.canonical) && attributes.size()==other.attributes.size()
				&& Objects.equals(attributes.fileKey(), other.attributes.fileKey())
				&& attributes.lastModifiedTime().equals(other.attributes.lastModifiedTime())
				&& attributes.creationTime().equals(other.attributes.creationTime());
		}
		final Path canonical;
		final BasicFileAttributes attributes;
	}

	private static final class Table{
		Table(String band_){band=band_;}
		int id(String label){
			require(label!=null && label.matches("[A-Z]+"), "Bad alphabet label: "+label);
			Integer id=map.get(label);
			if(id==null){id=labels.size(); labels.add(label); map.put(label, id);}
			return id;
		}
		void finish(){
			n=labels.size();
			require(n==FamilyDimerSubstitutionTable.DIMERS,
				"Rarity-aware sweep requires the canonical 400-symbol dimer dictionary in band "+band);
			verifyLabels(labels());
			separator=labels.get(0).length()==1 ? "" : ",";
			for(String label : labels){
				require(!separator.isEmpty() || label.length()==labels.get(0).length(),
					"Delimiter-free labels must have equal widths so partitions stay unambiguous");
			}
			odds=new double[n*n]; frequency=new double[n];
			final long[] endpoints=new long[n];
			final boolean[] endpointSeen=new boolean[n];
			final boolean[][] seen=new boolean[n][n];
			for(Cell cell : rows){
				require(cell.i<n && cell.j<n && !seen[cell.i][cell.j], "Duplicate symmetric cell");
				seen[cell.i][cell.j]=true;
				odds[cell.i*n+cell.j]=odds[cell.j*n+cell.i]=cell.odds;
				if(endpointSeen[cell.i]){require(endpoints[cell.i]==cell.endpointI, "Endpoint_i changed for "+labels.get(cell.i));}
				if(endpointSeen[cell.j]){require(endpoints[cell.j]==cell.endpointJ, "Endpoint_j changed for "+labels.get(cell.j));}
				require(cell.i!=cell.j || cell.endpointI==cell.endpointJ, "Diagonal endpoints must agree for "+labels.get(cell.i));
				endpoints[cell.i]=cell.endpointI; endpoints[cell.j]=cell.endpointJ;
				endpointSeen[cell.i]=endpointSeen[cell.j]=true;
			}
			long endpointTotal=0;
			for(int i=0; i<n; i++){endpointTotal=Math.addExact(endpointTotal, endpoints[i]);}
			require(endpointTotal>0, "Endpoint total must be positive");
			for(int i=0; i<n; i++){
				require(endpointSeen[i] && endpoints[i]>=0, "Every label needs a nonnegative endpoint total: "+i);
				frequency[i]=endpoints[i]/(double)endpointTotal;
				for(int j=i; j<n; j++){
					require(seen[i][j], "Missing upper-triangle cell "+i+","+j);
					require(Double.isFinite(odds[i*n+j]) && odds[i*n+j]>=0, "Bad odds "+i+","+j);
				}
			}
		}
		String[] labels(){return labels.toArray(new String[labels.size()]);}

		void verifyLabels(String[] array){
			require(Arrays.equals(array, FamilyDimerSubstitutionTable.dimerLabels()),
				"The 400-symbol source must use the canonical dimer dictionary order");
		}
		String separator(){return separator;}
		final String band;
		final ArrayList<String> labels=new ArrayList<String>();
		final HashMap<String,Integer> map=new HashMap<String,Integer>();
		final ArrayList<Cell> rows=new ArrayList<Cell>();
		int n;
		double[] odds, frequency;
		String separator;
	}

	private static final class Cell{
		Cell(int i_, int j_, long endpointI_, long endpointJ_, double odds_){i=i_; j=j_; endpointI=endpointI_; endpointJ=endpointJ_; odds=odds_;}
		final int i, j;
		final long endpointI, endpointJ;
		final double odds;
	}

	private static final class Group{
		Group(int[] members_, String text_, double frequency_){members=members_; text=text_; frequency=frequency_;}
		final int[] members;
		final String text;
		final double frequency;
	}

	private static final class Hierarchy{final ArrayList<Cut> cuts=new ArrayList<Cut>();}
	private static final class Cut{
		Cut(int symbols_, String[] groups_){symbols=symbols_; groups=groups_;}
		final int symbols;
		final String[] groups;
	}

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

	private static final String HEADER="band\tdimer_i\tdimer_j\tcount\tendpoint_i\tendpoint_j\tobserved_frequency\texpected_frequency\todds";
	private static final String PAIR_HEADER="view\tband\tquery_dimer\treference_dimer\tcount";
	private static final String BLEND_HEADER="dimer_i\tdimer_j\tconsensus_frequency\tpairwise_frequency\tblend50_frequency\tblend50_odds";
	private static final String BAND_BLEND_MAGIC="#bbgap_dimer_blend\t1";
	private static final String BAND_BLEND_HEADER="band\t"+BLEND_HEADER;
	private static final ArrayList<String> OPTIONS=new ArrayList<String>(Arrays.asList("in", "insha80", "out", "band", "lambdas", "targets", "verification", "format", "mergediagnostics"));
}
