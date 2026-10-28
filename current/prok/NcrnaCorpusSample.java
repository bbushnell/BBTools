package prok;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Locale;
import java.util.Random;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.IntHashMap;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Samples an annotated corpus without fitting models or changing annotations.
 * Probability arm: uniform loci from the declared proxy, excluding pilot source
 * groups. Challenge arm: round-robin strata, prioritizing proximity to the frozen
 * length-band boundaries within each stratum; distance ties are seeded. Source
 * probability=0 challenge=0 emits an unsampled census, including empty eligible
 * populations; it does not create a training/holdout assignment. Source
 * groups link taxid, versionless accession and identical compressed FASTA bytes;
 * they are independence blocks, not assertions that every linked genome is identical.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaCorpusSample {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args, NcrnaCorpusSample.class, false);
		String loci=null, genomes=null, hashes=null, pilot=null, out=null;
		int expected=-1, expectedGenomes=-1, low=-1, high=-1, probability=-1, challenge=-1;
		long seed=Long.MIN_VALUE;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0, "Expected flag=value: "+arg);
			final String key=arg.substring(0, eq).toLowerCase(Locale.ROOT), value=arg.substring(eq+1);
			if(key.equals("loci")){loci=value;}else if(key.equals("genomes")){genomes=value;}
			else if(key.equals("hashes")){hashes=value;}else if(key.equals("pilotgenomes")){pilot=value;}
			else if(key.equals("outprefix")){out=value;}else if(key.equals("expected")){expected=Parse.parseIntKMG(value);}
			else if(key.equals("expectedgenomes")){expectedGenomes=Parse.parseIntKMG(value);}
			else if(key.equals("minlen")){low=Parse.parseIntKMG(value);}else if(key.equals("maxlen")){high=Parse.parseIntKMG(value);}
			else if(key.equals("probability")){probability=Parse.parseIntKMG(value);}else if(key.equals("challenge")){challenge=Parse.parseIntKMG(value);}
			else if(key.equals("seed")){seed=Long.parseLong(value);}else{throw new IllegalArgumentException("Unknown parameter: "+key);}
		}
		require(loci!=null && genomes!=null && hashes!=null && pilot!=null && out!=null && expected>0 && expectedGenomes>0,
			"Bind the complete locus/genome/FASTA-hash ledgers and pilot roster with explicit denominators");
		require(low>=17 && high>=low && ((probability>0 && challenge>0) || (probability==0 && challenge==0)) && seed!=Long.MIN_VALUE,
			"Declare proxy bounds and seed; positive arms sample, both zero request a census without inventing eligible members");
		Shared.setThreads(1);
		final NcrnaCorpusSample tool=new NcrnaCorpusSample();
		tool.loadGenomes(genomes, hashes, pilot, expectedGenomes);
		tool.loadLoci(loci, expected, low, high);
		tool.select(probability, challenge, seed);
		tool.write(out, probability, challenge, seed, low, high);
		Shared.closeStream(pp.outstream);
	}

	/*--------------------------------------------------------------*/
	/*----------------         Input and groups     ----------------*/
	/*--------------------------------------------------------------*/

	private void loadGenomes(final String path, final String hashes, final String pilot, final int expected){
		assert(genomes.isEmpty()):"Build independence groups once over the full source universe";
		final ObjectMap<String,String> hashByPath=new ObjectMap<String,String>(String.class, String.class);
		ByteFile in=RrnaResourceIO.open(hashes);final LineParser1 p=new LineParser1('\t');
		try{for(byte[] row;(row=in.nextLine())!=null;){
			p.set(row);require(p.terms()==2, "FASTA provenance requires sha80 and absolute path");
			final String hash=p.parseString(0), file=p.parseString(1);
			require(hash.matches("[0-9a-f]{20}") && file.startsWith("/") && !hashByPath.contains(file), "Invalid or duplicate FASTA provenance: "+file);
			hashByPath.put(file, hash);
		}}finally{require(!in.close(), "FASTA hash ledger read failed");}
		in=RrnaResourceIO.open(path);
		try{
			RrnaResourceIO.header(in, GENOME_HEADER, path);
			for(byte[] row;(row=in.nextLine())!=null;){
				p.set(row);require(p.terms()==11, "Expected the annotated corpus11-field genome roster");
				final Genome g=new Genome(p.parseString(1), p.parseString(2), p.parseInt(3), p.parseString(4), p.parseString(5), p.parseString(7), p.parseInt(8));
				require(g.taxid>0 && !g.assembly.equals("NA") && g.assembly.lastIndexOf('.')>0 && g.loci>=0,
					"Source groups need explicit taxid/versioned assembly and nonnegative locus counts: "+g.name);
				require(!byGenome.contains(g.name), "Duplicate source genome: "+g.name);
				g.hash=hashByPath.get(g.fasta);require(g.loci==0 || g.hash!=null, "Positive source has no bound FASTA hash: "+g.name);
				byGenome.put(g.name, g);genomes.add(g);
			}
		}finally{require(!in.close(), "Genome roster read failed");}
		require(genomes.size()==expected, "Genome denominator changed: "+genomes.size()+" expected="+expected);
		Collections.sort(genomes, (a,b)->a.name.compareTo(b.name));
		parent=new int[genomes.size()];for(int i=0;i<parent.length;i++){parent[i]=i;genomes.get(i).index=i;}
		final IntHashMap taxids=new IntHashMap();
		final ObjectMap<String,Genome> assemblies=new ObjectMap<String,Genome>(String.class, Genome.class), bytes=new ObjectMap<String,Genome>(String.class, Genome.class);
		for(Genome g:genomes){
			if(taxids.contains(g.taxid)){union(g.index, taxids.get(g.taxid));}else{taxids.put(g.taxid, g.index);}
			link(assemblies, g.assembly.substring(0, g.assembly.lastIndexOf('.')), g);
			if(g.hash!=null){link(bytes, g.hash, g);}
		}
		in=RrnaResourceIO.open(pilot);int count=0;
		try{
			RrnaResourceIO.header(in, GENOME_HEADER, pilot);
			for(byte[] row;(row=in.nextLine())!=null;){
				p.set(row);require(p.terms()==11, "Pilot roster schema changed");
				final Genome g=byGenome.get(p.parseString(1));require(g!=null, "Pilot source absent from full roster");
				require(g.assembly.equals(p.parseString(2)) && g.taxid==p.parseInt(3) && g.fasta.equals(p.parseString(7)), "Pilot source identity disagrees with full roster");
				genomes.get(root(g.index)).pilot=true;count++;
			}
		}finally{require(!in.close(), "Pilot roster read failed");}
		require(count>0, "Probability validation cannot silently omit the pilot exclusion roster");
		for(Genome g:genomes){g.group=genomes.get(root(g.index));}
	}

	private void link(final ObjectMap<String,Genome> map, final String key, final Genome g){
		assert(key!=null && !key.isEmpty()):"Empty group keys would falsely unite unrelated sources";
		final Genome old=map.get(key);if(old==null){map.put(key, g);}else{union(g.index, old.index);}
	}
	private int root(int i){
		assert(i>=0 && i<parent.length):"Union-find indexes address the canonical sorted genome roster";
		while(parent[i]!=i){parent[i]=parent[parent[i]];i=parent[i];}return i;
	}
	private void union(final int a, final int b){
		final int x=root(a), y=root(b);parent[Math.max(x, y)]=Math.min(x, y);
		assert(root(a)==root(b)):"All source links must be transitive before pilot exclusion";
	}

	private void loadLoci(final String path, final int expected, final int low, final int high){
		assert(loci.isEmpty()):"Sample membership is chosen from one immutable source ledger";
		final ObjectSet<String> seen=new ObjectSet<String>(String.class);
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{
			RrnaResourceIO.header(in, NcrnaCmLocusAudit.HEADER, path);
			for(byte[] row;(row=in.nextLine())!=null;){
				p.set(row);require(p.terms()==21, "Locus ledger schema changed");
				final String id=p.parseString(0);require(seen.add(id), "Duplicate locus would bias sampling: "+id);
				final Genome g=byGenome.get(p.parseString(1));require(g!=null, "Locus source is not in the complete roster: "+id);
				require(g.taxid==p.parseInt(3) && g.assembly.equals(p.parseString(2)) && g.domain.equals(p.parseString(4)) && g.phylum.equals(p.parseString(5)), "Locus and genome identity disagree: "+id);
				final int len=p.parseInt(9), undefined=p.parseInt(17);require(len>0 && undefined>=0 && undefined<=len, "Locus lengths and ambiguity counts must be valid: "+id);
				final boolean edge=bool(p, 16), partial=bool(p, 12), pseudo=bool(p, 13), full=p.termEquals("ANNOTATION_FULL_COPY", 14);
				final boolean contiguous=p.parseString(8).indexOf(',')<0;
				require(!full || (!partial && !pseudo && contiguous), "Annotation-full state contradicts explicit source flags: "+id);
				final String bin=len<low?"SHORT":len>high?"LONG":"IN_RANGE";
				final String status=p.parseString(14);
				require(full || status.equals("SOURCE_PARTIAL") || status.equals("PSEUDOGENE") || status.equals("MULTIPART_UNRESOLVED"), "Unknown annotation eligibility state: "+status);
				final boolean eligible=full && !edge && undefined==0 && len>=low && len<=high;
				final String key=g.domain+"/"+g.phylum+"/"+status+"/"+bin+"/edge="+edge+"/ambiguous="+(undefined>0);
				Stratum stratum=byStratum.get(key);if(stratum==null){stratum=new Stratum(key);byStratum.put(key, stratum);strata.add(stratum);}
				final int distance=Math.min(Math.abs(len-low), Math.abs(len-high));
				final Locus l=new Locus(id, row, g, stratum, eligible, bin, distance);loci.add(l);stratum.population++;g.observed++;
				if(eligible){eligibleCount++;if(!g.group.pilot){probabilityPool.add(l);stratum.probabilityAvailable++;}}
			}
		}finally{require(!in.close(), "Locus ledger read failed");}
		require(loci.size()==expected, "Locus population changed: "+loci.size()+" expected="+expected);
		for(Genome g:genomes){require(g.observed==g.loci, "Genome/locus counts disagree: "+g.name);}
		Collections.sort(loci, (a,b)->a.id.compareTo(b.id));Collections.sort(probabilityPool, (a,b)->a.id.compareTo(b.id));
		Collections.sort(strata, (a,b)->a.key.compareTo(b.key));
	}

	/*--------------------------------------------------------------*/
	/*----------------         Select and write     ----------------*/
	/*--------------------------------------------------------------*/

	private void select(final int probability, final int challenge, final long seed){
		if(probability==0 && challenge==0){
			for(Locus l:loci){l.arm="CENSUS";l.stratum.remaining.add(l);}
			assert(!loci.isEmpty()):"A census preserves the validated complete source population even when no member is eligible";
			return;
		}
		require(probabilityPool.size()>=probability, "Not enough untouched proxy-eligible loci: "+probabilityPool.size()+" requested="+probability);
		final Random random=new Random(seed);Collections.shuffle(probabilityPool, random);
		for(int i=0;i<probability;i++){final Locus l=probabilityPool.get(i);l.arm="PROBABILITY";l.stratum.probabilitySelected++;}
		for(Locus l:loci){if(l.arm==null){l.stratum.remaining.add(l);}}
		int activeStrata=0;
		for(Stratum s:strata){
			Collections.shuffle(s.remaining, random);
			// Stable sorting preserves the seeded random order within each distance tie.
			Collections.sort(s.remaining, (a,b)->Integer.compare(a.distance, b.distance));
			if(!s.remaining.isEmpty()){activeStrata++;}
		}
		require(challenge>=activeStrata, "Challenge sample must cover every nonempty declared stratum: "+activeStrata+" requested="+challenge);
		int selected=0;
		while(selected<challenge){
			boolean progress=false;
			for(Stratum s:strata){
				if(selected==challenge){break;}
				if(s.challengeSelected>=s.remaining.size()){continue;}
				final Locus l=s.remaining.get(s.challengeSelected++);l.arm="CHALLENGE";selected++;progress=true;
			}
			require(progress, "Challenge population exhausted before the requested distinct sample size");
		}
		assert(selected==challenge):"Round-robin redistributes exhausted strata without silently changing the requested arm size";
	}

	private void write(final String prefix, final int probability, final int challenge, final long seed, final int low, final int high){
		final boolean census=probability==0 && challenge==0;
		RrnaResourceIO.fresh(prefix, ".ids.txt", ".eligible.ids.txt", ".loci.tsv", ".strata.tsv", ".groups.tsv", ".summary.tsv");
		final ByteStreamWriter ids=RrnaResourceIO.writer(prefix+".ids.txt"), rows=RrnaResourceIO.writer(prefix+".loci.tsv");
		final ByteStreamWriter eligible=RrnaResourceIO.writer(prefix+".eligible.ids.txt");
		final ByteStreamWriter groups=RrnaResourceIO.writer(prefix+".groups.tsv"), ss=RrnaResourceIO.writer(prefix+".strata.tsv"), summary=RrnaResourceIO.writer(prefix+".summary.tsv");
		final ByteBuilder b=new ByteBuilder();int n=0, writtenEligible=0;final ObjectSet<String> selectedGroups=new ObjectSet<String>(String.class);
		try{
			rows.println(NcrnaCmLocusAudit.HEADER+"\tsource_group\tpilot_exposed_group\tproxy_eligible\tlength_bin\tarm\tsplit\tstratum\tlength_boundary_distance");
			for(Locus l:loci){
				if(l.eligible){eligible.println(l.id);writtenEligible++;}
				if(l.arm!=null){
					require(!l.arm.equals("PROBABILITY") || (l.eligible && !l.genome.group.pilot), "Probability sampling must exclude all pilot-linked source groups");
					ids.println(l.id);selectedGroups.add(l.genome.group.name);n++;
					rows.print(b.clear().append(l.row).tab().append(l.genome.group.name).tab().append(l.genome.group.pilot).tab().append(l.eligible)
						.tab().append(l.bin).tab().append(l.arm).append("\tNOT_ASSIGNED\t").append(l.stratum.key).tab().append(l.distance).nl());
				}
			}
			require(writtenEligible==eligibleCount, "The full proxy population must be exported independently of sample and pilot exclusion");
			require(n==(census?loci.size():probability+challenge), "A census retains every source row; sample arms conserve their requested distinct counts");
			ss.println("stratum\tpopulation\tprobability_available\tprobability_selected\tchallenge_remaining\tchallenge_selected\tprobability_inclusion_numerator\tprobability_inclusion_denominator\tprobability_inclusion_applies_to");
			for(Stratum s:strata){ss.print(b.clear().append(s.key).tab().append(s.population).tab().append(s.probabilityAvailable).tab().append(s.probabilitySelected)
				.tab().append(s.remaining.size()).tab().append(s.challengeSelected).tab().append(s.probabilityAvailable==0?0:probability).tab().append(probabilityPool.size())
				.append("\tPROXY_ELIGIBLE_AND_NONPILOT_GROUP\n"));}
			groups.println("genome\tassembly\ttaxid\tfasta_sha80\tsource_group\tpilot_exposed_group");
			for(Genome g:genomes){groups.print(b.clear().append(g.name).tab().append(g.assembly).tab().append(g.taxid).tab().append(g.hash==null?"NA":g.hash).tab().append(g.group.name).tab().append(g.group.pilot).nl());}
			summary.println("population\tproxy_eligible\tprobability_available\tprobability_selected\tchallenge_selected\tselected_source_groups\tminlen\tmaxlen\tseed\tsplit_status\tchallenge_policy");
			summary.print(b.clear().append(loci.size()).tab().append(eligibleCount).tab().append(probabilityPool.size()).tab().append(probability).tab().append(challenge)
				.tab().append(selectedGroups.size()).tab().append(low).tab().append(high).tab().append(seed).append("\tNOT_ASSIGNED\t")
				.append(census?"CENSUS_NO_SAMPLING":"NEAREST_LENGTH_BOUNDARY_THEN_SEEDED_TIE").nl());
		}finally{RrnaResourceIO.finish(ids, eligible, rows, groups, ss, summary);}
		System.err.println("NCRNA_CORPUS_SAMPLE_PASS population="+loci.size()+" probability="+probability+" challenge="+challenge+"; no fitted resources or frozen train/holdout split");
	}

	private static boolean bool(final LineParser1 p, final int field){
		require(p.termEquals("true", field) || p.termEquals("false", field), "Source flags must be explicit true/false");return p.termEquals("true", field);
	}
	private static void require(final boolean ok, final String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final class Genome {
		Genome(String n, String a, int t, String d, String p, String f, int count){name=n;assembly=a;taxid=t;domain=d;phylum=p;fasta=f;loci=count;}
		final String name, assembly, domain, phylum, fasta;final int taxid, loci;
		String hash;int index, observed;boolean pilot;Genome group;
	}
	private static final class Locus {
		Locus(String i, byte[] r, Genome g, Stratum s, boolean e, String bin_, int distance_){id=i;row=r;genome=g;stratum=s;eligible=e;bin=bin_;distance=distance_;}
		final String id, bin;final byte[] row;final Genome genome;final Stratum stratum;final boolean eligible;final int distance;String arm;
	}
	private static final class Stratum {
		Stratum(String k){key=k;}final String key;final ArrayList<Locus> remaining=new ArrayList<Locus>();
		int population, probabilityAvailable, probabilitySelected, challengeSelected;
	}
	private final ObjectMap<String,Genome> byGenome=new ObjectMap<String,Genome>(String.class, Genome.class);
	private final ObjectMap<String,Stratum> byStratum=new ObjectMap<String,Stratum>(String.class, Stratum.class);
	private final ArrayList<Genome> genomes=new ArrayList<Genome>();
	private final ArrayList<Locus> loci=new ArrayList<Locus>(), probabilityPool=new ArrayList<Locus>();
	private final ArrayList<Stratum> strata=new ArrayList<Stratum>();
	private int[] parent;private int eligibleCount;
	public static final String GENOME_HEADER="index\tgenome\tassembly\ttaxid\tdomain\tphylum\tgff\tfasta\tloci\tgenome_bases\tcontigs";
}
