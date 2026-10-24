package prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Locale;
import java.util.Map;
import java.util.TreeMap;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;
import prok.NcrnaGffFamilyGrader.Locus;
import prok.NcrnaGffFamilyGrader.MatchResult;
import prok.NcrnaGffFamilyGrader.MatchedPair;

/** Annotation-relative grading; no independent biological truth claim.
 * Reuses the reviewed maximum-cardinality matcher separately per contig/family.
 * Declared circular source features retain their unrolled coordinates. Circular
 * matching counts physical overlap while retaining one truth index; endpoints
 * use nearest signed modular offsets and frames use the lifted5prime anchor.
 * Linear contigs retain the original metric. No annotation is clipped or dropped.
 * The run producer must verify the declared input and GFF hashes before invocation.
 * @author Ganyu */
public final class PanelAnnotationGrade {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args, PanelAnnotationGrade.class, false);
		String roster=null, references=null, arms=null, out=null, bacteria=null, archaea=null, baseline="default";
		int expected=203, expectedRoster=223, expectedBacteria=100, expectedArms=4;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0, "Expected flag=value: "+arg);
			final String key=arg.substring(0, eq).toLowerCase(Locale.ROOT), value=arg.substring(eq+1);
			if(key.equals("roster")){roster=value;}else if(key.equals("references")){references=value;}
			else if(key.equals("arms")){arms=value;}else if(key.equals("outprefix")){out=value;}
			else if(key.equals("bacteriaroot")){bacteria=value;}else if(key.equals("archaearoot")){archaea=value;}
			else if(key.equals("baseline")){baseline=value;}else if(key.equals("expected")){expected=Integer.parseInt(value);}
			else if(key.equals("expectedroster")){expectedRoster=Integer.parseInt(value);}
			else if(key.equals("expectedbacteria")){expectedBacteria=Integer.parseInt(value);}
			else if(key.equals("expectedarms")){expectedArms=Integer.parseInt(value);}
			else{throw new IllegalArgumentException("Unknown flag: "+key);}
		}
		require(roster!=null && references!=null && arms!=null && out!=null && bacteria!=null && archaea!=null,
			"Bind roster, reference manifest, arm manifest, both explicit source roots and fresh output");
		require(expectedArms>0,"Expected arm count must be positive");
		require(expected>0 && expectedRoster>=expected && expectedBacteria>=0 && expectedBacteria<=expected,
			"Expected counts must describe the declared reference subset");
		require(bacteria.endsWith("/") && archaea.endsWith("/") && !bacteria.startsWith(archaea) && !archaea.startsWith(bacteria),
			"Domain membership uses disjoint explicit source-directory prefixes, never taxid/name inference");
		Shared.setThreads(1);ReadWrite.setZipThreads(1);ReadWrite.ZIPLEVEL=2;
		ReadWrite.USE_BGZF=true;ReadWrite.ALLOW_NATIVE_BGZF=true;ReadWrite.PREFER_BGZIP=true;
		final PanelAnnotationGrade grader=new PanelAnnotationGrade();
		grader.loadRoster(roster, expectedRoster);
		grader.loadReferences(references, bacteria, archaea, expected, expectedBacteria);
		grader.loadContigs();grader.loadArms(arms, baseline, expectedArms);grader.run(out);
		pp.outstream.println("PANEL_ANNOTATION_GRADE_COMPLETE reference_genomes="+expected+" arms="+grader.arms.size()+"; annotation_relative=true; coordinate_metric=EXPLICIT_CIRCULAR_ELSE_LINEAR; origin_spanning_source_rows="+grader.originRows);
		Shared.closeStream(pp.outstream);
	}

	void loadRoster(String path, int expected){
		assert(genomes.isEmpty()):"Load one source roster before annotation or caller evidence";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{
			RrnaResourceIO.header(in, "clade\tassembly\tfasta\tfasta_sha80", path);
			for(byte[] line;(line=in.nextLine())!=null;){
				p.set(line);require(p.terms()==4, "Roster must have four fields");
				final Genome g=new Genome(p.parseString(0), p.parseString(1), p.parseString(2), p.parseString(3));
				require(!g.id.isEmpty() && SeedOffsetTable.validSha80(g.fastaHash) && byGenome.put(g.id, g)==null,
					"Roster identifiers must be unique with valid source hashes");genomes.add(g);
			}
		}finally{require(!in.close(), "Roster read failed");}
		require(genomes.size()==expected, "Declared caller roster count differs: "+genomes.size());
	}
	void loadReferences(String path, String bacteria, String archaea, int expected, int expectedBacteria){
		assert(!genomes.isEmpty() && referenceGenomes.isEmpty()):"References must bind a subset of the caller roster";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');int bacterial=0;
		try{
			RrnaResourceIO.header(in, "index\tgenome\tfasta\tfasta_sha80\tgff\tgff_sha80", path);
			for(byte[] line;(line=in.nextLine())!=null;){
				p.set(line);require(p.terms()==6 && p.parseInt(0)==referenceGenomes.size()+1, "Reference indices must be sequential");
				final Genome g=byGenome.get(p.parseString(1));
				require(g!=null && g.referenceIndex<0 && p.termEquals(g.fasta, 2) && p.termEquals(g.fastaHash, 3),
					"Reference assembly/path/hash must exactly match caller roster");
				require(g.clade.equals("prok"), "Reference membership must be declared prok in source roster");
				g.domain=g.fasta.startsWith(bacteria)?BACTERIA:g.fasta.startsWith(archaea)?ARCHAEA:-1;
				require(g.domain>=0, "Reference source is outside both declared domain roots: "+g.fasta);
				if(g.domain==BACTERIA){bacterial++;}
				g.gff=p.parseString(4);g.gffHash=p.parseString(5);
				require(SeedOffsetTable.validSha80(g.gffHash), "Reference GFF hash is required");
				g.referenceIndex=referenceGenomes.size();referenceGenomes.add(g);
			}
		}finally{require(!in.close(), "Reference manifest read failed");}
		require(referenceGenomes.size()==expected && bacterial==expectedBacteria, "Reference/domain membership changed");
	}
	void loadContigs(){
		assert(contigs.size()==0):"The actual caller FASTA set establishes unique literal seqid identity once";
		for(Genome g:genomes){
			RrnaResourceIO.stream(g.fasta, read->{
				final String seqid=RrnaResourceIO.first(read.id);require(read.length()>0, "Empty caller contig: "+seqid);
				final Contig c=new Contig(g, seqid, read.length());
				require(contigs.put(seqid, c)==null, "Repeated native seqid across caller inputs needs explicit disambiguation: "+seqid);
				g.contigs.add(c);
			});
			require(!g.contigs.isEmpty(), "Empty roster genome: "+g.id);
		}
	}
	void loadArms(String path, String baseline){
		loadArms(path,baseline,4);
	}
	void loadArms(String path, String baseline, int expected){
		assert(arms.isEmpty()):"Load one configuration manifest";
		require(expected>0,"Expected arm count must be positive");
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		final ObjectSet<String> names=new ObjectSet<String>(String.class);Arm control=null;
		try{
			RrnaResourceIO.header(in, "arm\tgff\tgff_sha80", path);
			for(byte[] line;(line=in.nextLine())!=null;){
				p.set(line);require(p.terms()==3, "Arm manifest requires three fields");
				final Arm a=new Arm(p.parseString(0), p.parseString(1), p.parseString(2));
				require(!a.name.isEmpty() && names.add(a.name) && SeedOffsetTable.validSha80(a.hash), "Unique named arms and GFF hashes required");
				arms.add(a);if(a.name.equals(baseline)){control=a;}
			}
		}finally{require(!in.close(), "Arm manifest read failed");}
		require(arms.size()==expected && control!=null, "Expected "+expected+" arms including explicit baseline; observed "+arms.size());
		arms.remove(control);arms.add(0, control); // Baseline coordinate multisets are retained for all later comparisons.
	}

	void readGff(String path, Genome truthGenome, Arm arm, ByteStreamWriter counts, ByteStreamWriter ungradedRows, ByteStreamWriter circularRows){
		assert((truthGenome==null)!=(arm==null)):"A GFF is either source annotation or one caller arm";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		final TreeMap<String,long[][]> typeCounts=new TreeMap<String,long[][]>();long lineNumber=0, features=0, counted=0, ungraded=0, expectedUngraded=0;
		final ByteBuilder b=new ByteBuilder();
		String lastSeq=null;Contig current=null;
		try{for(byte[] line;(line=in.nextLine())!=null;){
			lineNumber++;if(line.length==0){continue;}
			if(line[0]=='#'){if(shared.Tools.startsWith(line, "##FASTA")){break;}continue;}
			p.set(line);require(p.terms()==9, "Expected GFF9 at "+path+":"+lineNumber);
			if(lastSeq==null || !p.termEquals(lastSeq, 0)){lastSeq=p.parseString(0);current=contigs.get(lastSeq);}
			require(current!=null && (truthGenome==null || current.genome==truthGenome), "GFF seqid is not in its bound source genome: "+lastSeq);
			final int start=p.parseInt(3), stop=p.parseInt(4);
			final Genome g=current.genome;final String type=featureType(p);
			validateInterval(current, start, stop, truthGenome!=null, path+":"+lineNumber);
			if(truthGenome!=null && type.equals("region")){
				final String circular=RrnaResourceIO.optionalField(p.parseString(8), "Is_circular");
				if(circular!=null){
					require(start==1 && stop==current.length && (circular.equals("true") || circular.equals("false")),
						"Circular topology must be explicit on a full actual-length region: "+path+":"+lineNumber);
					final boolean value=circular.equals("true");
					require(!current.topologyKnown || current.circular==value, "Conflicting topology for "+current.id);
					current.topologyKnown=true;current.circular=value;
				}
			}
			final int family=family(type, p, stop-start+1, truthGenome!=null);
			final int scope=g.referenceIndex<0?2:g.domain;
			if(stop>current.length){
				current.wrappedSource=true;originRows++;
				b.clear().append(g.id).tab().append(DOMAIN[g.domain]).tab().append(current.length).tab().append(lineNumber)
					.tab().append(family<0?"UNGRADED":FAMILY[family]).tab().append(line).nl();circularRows.print(b);
			}
			String countKey=type;
			if(family<0 && (type.equals("RNA") || type.equals("rRNA"))){
				final String attrs=p.parseString(8), tag=colonField(attrs,"family"), model=colonField(attrs,"model");
				countKey=type+'\t'+(tag==null?"NA":tag)+'\t'+(model==null?"NA":model);
			}
			long[][] n=typeCounts.get(countKey);if(n==null){n=new long[3][FAMILY.length+1];typeCounts.put(countKey, n);}n[scope][family<0?FAMILY.length:family]++;features++;
			if(arm!=null && (g.referenceIndex<0 || family<0) && (type.equals("RNA") || type.equals("rRNA"))){
				b.clear().append(arm.name).tab().append(g.id).tab().append(scope==2?"OUTSIDE_REFERENCE":DOMAIN[scope]).tab().append(lineNumber).tab().append(line).nl();
				ungradedRows.print(b);ungraded++;
			}
			if(g.referenceIndex<0 || family<0){continue;}
			final char strand=p.termEquals("+",6)?'+':p.termEquals("-",6)?'-':'?';
			require(strand!='?', "Graded features need known strand at "+path+":"+lineNumber);
			final int phase=family!=CDS || p.termEquals(".",7)?-1:p.termEquals("0",7)?0:p.termEquals("1",7)?1:p.termEquals("2",7)?2:-2;
			require(phase>=-1, "CDS phase must be0/1/2 or missing");
			final Feature f=new Feature(new Locus(current.id, FAMILY[family], start, stop, strand), phase, lineNumber);
			(truthGenome==null?current.calls[family]:current.truth[family]).add(f);
		}}finally{require(!in.close(), "GFF read failed: "+path);}
		//Topology can occur after a feature in a valid source GFF; validate before any matching.
		if(truthGenome!=null){for(Contig c:truthGenome.contigs){
			require(!c.wrappedSource || c.topologyKnown && c.circular,
				"Origin-spanning annotation lacks an explicit full-length circular region: "+path+" "+c.id);
		}}
		for(Map.Entry<String,long[][]> e:typeCounts.entrySet()){
			final String key=e.getKey();final boolean rna=key.equals("RNA") || key.equals("rRNA") || key.startsWith("RNA\t") || key.startsWith("rRNA\t");
			for(int scope=0;scope<3;scope++){for(int family=0;family<=FAMILY.length;family++){
				final long n=e.getValue()[scope][family];if(n==0){continue;}counted+=n;
				if(arm!=null && rna && (scope==2 || family==FAMILY.length)){expectedUngraded+=n;}
				b.clear().append(arm==null?"REFERENCE":arm.name).tab().append(truthGenome==null?"CALLER_COHORT":truthGenome.id)
					.tab().append(scope==2?"OUTSIDE_REFERENCE":DOMAIN[scope]).tab().append(e.getKey());
				if(e.getKey().indexOf('\t')<0){b.append("\tNA\tNA");}
				b.tab().append(family==FAMILY.length?"UNGRADED":FAMILY[family]).tab().append(n).nl();counts.print(b);
			}}
		}
		require(counted==features, "Every GFF feature must appear in declared-scope type counts");
		require(ungraded==expectedUngraded, "Every ungraded RNA/rRNA call must retain its coordinates and original nine fields");
	}
	public static void validateInterval(Contig c, int start, int stop, boolean source, String location){
		require(start>=1 && start<=c.length && stop>=start && (long)stop-start+1<=c.length && (source || stop<=c.length),
			"GFF interval is outside actual contig or spans multiple turns: "+location);
	}
	static int family(String type, LineParser1 p, int length, boolean reference){
		assert(type!=null):"Attribution is by explicit GFF type and established rRNA parser";
		if(type.equals("CDS")){return CDS;}if(type.equals("tRNA")){return TRNA;}
		if(!type.equals("rRNA")){return -1;}
		final String attrs=p.parseString(8), tag=colonField(attrs,"family");
		if(!reference && tag!=null){
			if(tag.equals("prok5S")){return FIVE;}if(tag.equals("16S")){return SSU;}if(tag.equals("23S")){return LSU;}
			return -1; // An explicitly tagged euk5S is never a prok5S call, even if its display label is5S.
		}
		if(!reference){final String model=colonField(attrs,"model");if(model!=null && model.startsWith("euk5S_")){return -1;}}
		final String name=LegacyRrnaGffGrader.extractLegacyFamily(attrs, length);
		if("5S".equals(name)){return FIVE;}if("16S".equals(name)){return SSU;}if("23S".equals(name)){return LSU;}return -1;
	}
	/** Exact BBTools comma/semicolon-delimited colon fields, not substring attribution. */
	static String colonField(String attrs,String key){
		assert(attrs!=null && key!=null):"Only explicit attribute fields may override displayed rRNA names";
		final String prefix=key+':';int from=0;
		while(from<attrs.length()){
			final int i=attrs.indexOf(prefix,from);if(i<0){return null;}
			if(i==0 || attrs.charAt(i-1)==',' || attrs.charAt(i-1)==';'){
				int end=i+prefix.length();while(end<attrs.length() && attrs.charAt(end)!=',' && attrs.charAt(end)!=';'){end++;}
				final String value=attrs.substring(i+prefix.length(),end).trim();require(!value.isEmpty(),"Empty explicit "+key+" field");return value;
			}from=i+prefix.length();
		}return null;
	}
	static String featureType(LineParser1 p){
		for(String type:COMMON_TYPES){if(p.termEquals(type,2)){return type;}}
		return p.parseString(2); // Only uncommon types allocate a retained-count-map key candidate.
	}

	void run(String prefix){
		assert(!arms.isEmpty() && !referenceGenomes.isEmpty()):"Bound callers must be graded against one fixed source subset";
		RrnaResourceIO.fresh(prefix, ".summary.tsv", ".deltas.tsv", ".loci.tsv.gz", ".type_counts.tsv", ".ungraded_rna.tsv.gz", ".circular_source.tsv.gz");
		final ByteStreamWriter summary=RrnaResourceIO.writer(prefix+".summary.tsv"), delta=RrnaResourceIO.writer(prefix+".deltas.tsv"),
			loci=RrnaResourceIO.writer(prefix+".loci.tsv.gz"), counts=RrnaResourceIO.writer(prefix+".type_counts.tsv"), ungraded=RrnaResourceIO.writer(prefix+".ungraded_rna.tsv.gz"),
			circular=RrnaResourceIO.writer(prefix+".circular_source.tsv.gz");
		final long[][][][] all=new long[arms.size()][referenceGenomes.size()][FAMILY.length][NSTAT];
		try{
			summary.println("arm\tscope\tfamily\t"+STAT_HEADER+"\tprecision\trecall\texact5_fraction\texact3_fraction\tboth_fraction\tframe_agreement_fraction\ttruth_source");
			delta.println("arm\tbaseline\tscope\tfamily\t"+STAT_HEADER+"\tprecision_delta\trecall_delta\tboth_fraction_delta");
			loci.println("arm\tgenome\tfamily\trole\tseqid\tstart\tstop\tstrand\tphase\tsource_line\toutcome\tpartner_source_line\tpartner_start\tpartner_stop\tpartner_phase\td5\td3\tcircular_length_modular_ends_frame");
			counts.println("arm\tgenome_or_cohort\tscope\toriginal_type\texplicit_family_tag\tmodel_tag\trecognized_grading_family\tfeature_rows");
			ungraded.println("arm\tgenome\tscope\tsource_line\toriginal_seqid\toriginal_source\toriginal_type\toriginal_start\toriginal_stop\toriginal_score\toriginal_strand\toriginal_phase\toriginal_attributes");
			circular.println("genome\tscope\tcontig_length\tsource_line\tgrading_family\toriginal_seqid\toriginal_source\toriginal_type\toriginal_start\toriginal_stop\toriginal_score\toriginal_strand\toriginal_phase\toriginal_attributes");
			for(Genome g:referenceGenomes){readGff(g.gff, g, null, counts, ungraded, circular);}
			for(int a=0;a<arms.size();a++){
				final Arm arm=arms.get(a);readGff(arm.gff, null, arm, counts, ungraded, circular);
				for(Genome g:referenceGenomes){
					for(Contig c:g.contigs){for(int f=0;f<FAMILY.length;f++){
						final int circle=c.circular?c.length:0;
						final MatchResult m=match(c.calls[f], c.truth[f], circle);
						final long[] s=statistics(c.calls[f], c.truth[f], m, f==CDS, circle);
						if(a>0){changes(c.calls[f], c.baseline[f], s);}
						add(all[a][g.referenceIndex][f], s);writeLoci(loci, arm.name, g.id, f, c.calls[f], c.truth[f], m, circle);
					}}
					for(int f=0;f<FAMILY.length;f++){
						writeSummary(summary, arm.name, g.id, f, all[a][g.referenceIndex][f]);
						writeDelta(delta, arm.name, arms.get(0).name, g.id, f, all[a][g.referenceIndex][f], all[0][g.referenceIndex][f]);
					}
					for(Contig c:g.contigs){if(a==0){c.baseline=c.calls;}c.calls=buckets();}
				}
				for(int domain=-1;domain<=ARCHAEA;domain++){
					final String scope=domain<0?"ALL_REFERENCE":DOMAIN[domain];
					for(int f=0;f<FAMILY.length;f++){
						final long[] s=new long[NSTAT], b=new long[NSTAT];
						for(Genome g:referenceGenomes){if(domain<0 || g.domain==domain){add(s, all[a][g.referenceIndex][f]);add(b, all[0][g.referenceIndex][f]);}}
						writeSummary(summary, arm.name, scope, f, s);writeDelta(delta, arm.name, arms.get(0).name, scope, f, s, b);
					}
				}
			}
		}finally{RrnaResourceIO.finish(summary, delta, loci, counts, ungraded, circular);}
	}
	public static MatchResult match(ArrayList<Feature> calls, ArrayList<Feature> truth){
		return match(calls, truth, 0);
	}
	public static MatchResult match(ArrayList<Feature> calls, ArrayList<Feature> truth, int circle){
		assert(calls!=null && truth!=null):"Do not pool contigs or families into a quadratic whole-cohort graph";
		final ArrayList<Locus> c=new ArrayList<Locus>(calls.size()), t=new ArrayList<Locus>(truth.size());
		for(Feature f:calls){c.add(f.locus);}for(Feature f:truth){t.add(f.locus);}
		return NcrnaGffFamilyGrader.matchOneToOne(c, t, .5f, circle);
	}
	public static long[] statistics(ArrayList<Feature> calls, ArrayList<Feature> truth, MatchResult m, boolean cds){
		return statistics(calls, truth, m, cds, 0);
	}
	public static long[] statistics(ArrayList<Feature> calls, ArrayList<Feature> truth, MatchResult m, boolean cds, int circle){
		assert(m.callMatchedTruthIdx.length==calls.size() && m.truthMatchedCallIdx.length==truth.size()):"Pair indices belong to these exact ordered lists";
		final long[] s=new long[NSTAT];s[0]=m.matched;s[1]=calls.size()-m.matched;s[2]=truth.size()-m.matched;
		s[6]=duplicates(calls);s[7]=duplicates(truth);
		for(int i=0;i<calls.size();i++){
			final int partner=m.callMatchedTruthIdx[i];if(partner<0){continue;}
			final Feature c=calls.get(i), t=truth.get(partner);final MatchedPair p=new MatchedPair(c.locus.family, t.locus, c.locus);
			final int d5=endpointDelta(p.d5, circle), d3=endpointDelta(p.d3, circle);
			if(d5==0){s[3]++;}if(d3==0){s[4]++;}if(d5==0 && d3==0){s[5]++;}
			if(cds){if(c.phase<0 || t.phase<0){s[10]++;}
				else if(circle==0 ? frame(c)==frame(t) : Math.floorMod((long)d5+c.phase-t.phase,3L)==0){s[8]++;}else{s[9]++;}}
		}
		require(!cds || s[8]+s[9]+s[10]==s[0], "Every CDS TP needs an explicit frame outcome");return s;
	}
	static int frame(Feature f){
		assert(f.phase>=0 && f.phase<=2):"Missing GFF phases must not become invented frame agreement";
		final long anchor=f.locus.strand=='+'?(long)f.locus.start+f.phase:(long)f.locus.stop-f.phase;
		return (int)Math.floorMod(anchor, 3L);
	}
	/** Nearest signed sense-coordinate offset, with negative half-turn ties.
	 * Circular CDS frames lift the call's5prime anchor into this reference-relative
	 * frame before applying GFF phases; a circumference need not be divisible by3. */
	public static int endpointDelta(int raw, int circle){
		assert(circle>=0):"Zero means the unchanged linear endpoint contract";
		return circle==0?raw:(int)(Math.floorMod((long)raw+circle/2,(long)circle)-circle/2);
	}
	static Feature[] sorted(ArrayList<Feature> list){
		final Feature[] a=list.toArray(new Feature[list.size()]);Arrays.sort(a, COORDINATES);return a;
	}
	static long duplicates(ArrayList<Feature> list){
		final Feature[] a=sorted(list);long n=0;for(int i=1;i<a.length;i++){if(COORDINATES.compare(a[i-1], a[i])==0){n++;}}return n;
	}
	static void changes(ArrayList<Feature> calls, ArrayList<Feature> baseline, long[] s){
		assert(s.length==NSTAT):"Coordinate deltas retain multiplicity, separately from source-annotation matching";
		final Feature[] c=sorted(calls), b=sorted(baseline);int i=0,j=0;
		while(i<c.length && j<b.length){final int order=COORDINATES.compare(c[i], b[j]);if(order==0){i++;j++;}else if(order<0){s[11]++;i++;}else{s[12]++;j++;}}
		s[11]+=c.length-i;s[12]+=b.length-j;
		require(c.length-b.length==s[11]-s[12], "Coordinate multiset delta must conserve caller count");
	}
	static void writeSummary(ByteStreamWriter out, String arm, String scope, int family, long[] s){
		final ByteBuilder b=new ByteBuilder().append(arm).tab().append(scope).tab().append(FAMILY[family]);
		for(long n:s){b.tab().append(n);}fraction(b,s[0],s[0]+s[1]);fraction(b,s[0],s[0]+s[2]);
		fraction(b,s[3],s[0]);fraction(b,s[4],s[0]);fraction(b,s[5],s[0]);fraction(b,s[8],s[8]+s[9]);out.print(b.append("\tSOURCE_GFF_FEATURE_ROWS").nl());
	}
	static void writeDelta(ByteStreamWriter out, String arm, String baseline, String scope, int family, long[] s, long[] control){
		final ByteBuilder b=new ByteBuilder().append(arm).tab().append(baseline).tab().append(scope).tab().append(FAMILY[family]);
		for(int i=0;i<NSTAT;i++){b.tab().append(s[i]-control[i]);}
		difference(b,s[0],s[0]+s[1],control[0],control[0]+control[1]);difference(b,s[0],s[0]+s[2],control[0],control[0]+control[2]);
		difference(b,s[5],s[0],control[5],control[0]);out.print(b.nl());
	}
	public static void writeLoci(ByteStreamWriter out, String arm, String genome, int family, ArrayList<Feature> calls, ArrayList<Feature> truth, MatchResult m, int circle){
		final ByteBuilder b=new ByteBuilder();
		for(int role=0;role<2;role++){
			final ArrayList<Feature> list=role==0?calls:truth, other=role==0?truth:calls;
			for(int i=0;i<list.size();i++){
				final Feature f=list.get(i);final Locus x=f.locus;final int partner=role==0?m.callMatchedTruthIdx[i]:m.truthMatchedCallIdx[i];
				b.clear().append(arm).tab().append(genome).tab().append(FAMILY[family]).tab().append(role==0?"CALL":"REFERENCE").tab()
					.append(x.seqid).tab().append(x.start).tab().append(x.stop).tab().append(x.strand).tab().append(f.phase).tab().append(f.line).tab()
					.append(partner>=0?"TP":role==0?"FP":"FN");
				if(partner<0){b.append("\tNA\tNA\tNA\tNA\tNA\tNA");}
				else{
					final Feature y=other.get(partner);final MatchedPair p=new MatchedPair(FAMILY[family],role==0?y.locus:x,role==0?x:y.locus);
						b.tab().append(y.line).tab().append(y.locus.start).tab().append(y.locus.stop).tab().append(y.phase)
							.tab().append(endpointDelta(p.d5,circle)).tab().append(endpointDelta(p.d3,circle));
					}out.print(b.tab().append(circle).nl());
			}
		}
	}
	static void add(long[] a,long[] b){assert(a.length==b.length):"Same metric schema required";for(int i=0;i<a.length;i++){a[i]+=b[i];}}
	static void fraction(ByteBuilder b,long n,long d){b.tab();if(d==0){b.append("NA");}else{b.append(n/(double)d,8);}}
	static void difference(ByteBuilder b,long n,long d,long bn,long bd){b.tab();if(d==0 || bd==0){b.append("NA");}else{b.append(n/(double)d-bn/(double)bd,8);}}
	public static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	@SuppressWarnings("unchecked") static ArrayList<Feature>[] buckets(){final ArrayList<Feature>[] a=new ArrayList[FAMILY.length];for(int i=0;i<a.length;i++){a[i]=new ArrayList<Feature>();}return a;}
	static final Comparator<Feature> COORDINATES=(a,b)->{int c=Character.compare(a.locus.strand,b.locus.strand);if(c==0){c=Integer.compare(a.locus.start,b.locus.start);}return c==0?Integer.compare(a.locus.stop,b.locus.stop):c;};
	public static final int CDS=0,TRNA=1,FIVE=2,SSU=3,LSU=4,BACTERIA=0,ARCHAEA=1,NSTAT=13;
	static final String[] FAMILY={"CDS","tRNA","5S","16S","23S"},DOMAIN={"BACTERIA","ARCHAEA"};
	static final String[] COMMON_TYPES={"CDS","tRNA","rRNA","RNA","gene","region","exon","ncRNA","misc_feature","repeat_region","tmRNA"};
	static final String STAT_HEADER="TP\tFP\tFN\texact5\texact3\tbothExact\tduplicateCalls\tduplicateReference\tCDS_frame_agree\tCDS_frame_disagree\tCDS_frame_unassessable\tcoordinate_added_vs_baseline\tcoordinate_removed_vs_baseline";
	public static final class Feature{public Feature(Locus l,int p,long n){locus=l;phase=p;line=n;}public final Locus locus;final int phase;final long line;}
	public static final class Contig{public Contig(Genome g,String i,int n){genome=g;id=i;length=n;}final Genome genome;final String id;final int length;boolean topologyKnown,circular,wrappedSource;final ArrayList<Feature>[] truth=buckets();ArrayList<Feature>[] calls=buckets(),baseline;}
	static final class Genome{Genome(String c,String i,String f,String h){clade=c;id=i;fasta=f;fastaHash=h;}final String clade,id,fasta,fastaHash;String gff,gffHash;int referenceIndex=-1,domain=-1;final ArrayList<Contig> contigs=new ArrayList<Contig>();}
	static final class Arm{Arm(String n,String g,String h){name=n;gff=g;hash=h;}final String name,gff,hash;}
	final ArrayList<Genome> genomes=new ArrayList<Genome>(),referenceGenomes=new ArrayList<Genome>();
	final ArrayList<Arm> arms=new ArrayList<Arm>();
	final ObjectMap<String,Genome> byGenome=new ObjectMap<String,Genome>(String.class,Genome.class);
	final ObjectMap<String,Contig> contigs=new ObjectMap<String,Contig>(String.class,Contig.class);
	long originRows=0;
}
