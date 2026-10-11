package prot;

import java.nio.charset.StandardCharsets;
import java.util.HashMap;
import java.util.TreeMap;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;
import structures.IntHashMap;

/** Complete-genome family presence and copy statistics for v1.4 inputs.
 * Populations use the same accepted QuickClade labels as expected-copy tables.
 * Stored integer numerators retain carrier-based and whole-population meanings.
 * @author Brian Bushnell, Yoimiya
 */
public final class MagQCRequirednessTable {

	/** Builds one immutable resource from explicitly pinned, existing corpus inputs. */
	public static void main(String[] args){
		final PreParser pp=new PreParser(args, MagQCRequirednessTable.class, false);
		try{
			final HashMap<String,String> options=new HashMap<String,String>();
			final Parser common=new Parser();
			for(String arg:pp.args){
				final int eq=arg.indexOf('=');
				require(eq>0, "Expected flag=value");
				final String key=arg.substring(0, eq).toLowerCase(java.util.Locale.ROOT), value=arg.substring(eq+1);
				if(key.equals("out") || key.equals("expectedunits") || key.equals("mode") || key.equals("count")
					|| key.equals("requiredness") || key.equals("requirednesssha80") || isInput(key)){
					require(!value.isEmpty() && options.put(key, value)==null, "Empty or duplicate table option: "+key);
				}else{require(common.parse(arg, key, value), "Unknown table option: "+key);}
			}
			if("markers".equals(options.get("mode"))){writeMarkers(options); return;}
			require(!options.containsKey("mode") || "build".equals(options.get("mode")), "Unknown requiredness mode");
			for(String name:INPUTS){
				require(DigestSuffix.file(required(options, name)).equals(DigestSuffix.requireSuffix(required(options, name+"sha80"), name)),
					"Requiredness input pin mismatch: "+name);
			}
			build(options);
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Publishes the requested universal single-copy ranks, retaining counts and denominator for audit. */
	private static void writeMarkers(final HashMap<String,String> options){
		final String family=required(options, "families"), familyPin=required(options, "familiessha80");
		require(DigestSuffix.file(family).equals(DigestSuffix.requireSuffix(familyPin, "families")), "Family input pin differs");
		final String source=required(options, "requiredness"), sourcePin=required(options, "requirednesssha80");
		final MagQCRequirednessTable table=load(source, sourcePin, familyCount(family), familyPin);
		final int[] ranks=table.topSingletonRanks(Integer.parseInt(required(options, "count")));
		final ByteBuilder out=new ByteBuilder("#schema_version\tmagqc_size_markers_v1\n#requiredness_sha80\t")
			.append(sourcePin).append("\n#family_sha80\t").append(familyPin)
			.append("\n#selection\tsingle_copy_prevalence_descending_then_family_rank_ascending\n")
			.append("#columns\torder\tfamily_rank\tsingle_copy_genomes\tcomplete_genomes\n");
		for(int i=0; i<ranks.length; i++){
			out.append(i).tab().append(ranks[i]).tab().append(table.global.singletons[ranks[i]]).tab().append(table.global.genomes).nl();
		}
		try(OutputFileSetTransaction tx=new OutputFileSetTransaction(new String[]{required(options, "out")}, new String[]{source,family}, false)){
			final ByteStreamWriter writer=new ByteStreamWriter(tx.temporaryPath(0), false, false, false);
			writer.start(); writer.print(out); require(!writer.poisonAndWait(), "Size marker resource write failed"); tx.commit();
		}
		System.err.println("SIZE_MARKERS_PASS count="+ranks.length+" genomes="+table.global.genomes);
	}

	/** Selects by singleton count since every family has the same complete-corpus denominator. */
	int[] topSingletonRanks(final int count){
		require(count>0 && count<=families, "Marker count must be positive and no larger than the family vocabulary");
		final int[] ranks=new int[count]; final boolean[] chosen=new boolean[families];
		for(int order=0; order<count; order++){
			int best=-1;
			for(int rank=0; rank<families; rank++){
				if(!chosen[rank] && (best<0 || global.singletons[rank]>global.singletons[best])){best=rank;}
			}
			assert(best>=0) : "Bounded top-K selection must find an unselected family at every output order";
			ranks[order]=best; chosen[best]=true;
		}
		return ranks;
	}

	/** Joins only the clean roster, preserving every eligible zero-family genome in denominators. */
	private static void build(HashMap<String,String> o){
		assert(o!=null) : "The CLI must validate the pinned resource arguments before building";
		final int families=familyCount(required(o, "families"));
		final HashMap<String,Genome> byUnit=new HashMap<String,Genome>();
		final IntHashMap byTid=new IntHashMap();
		final java.util.ArrayList<Genome> genomes=new java.util.ArrayList<Genome>();
		final LineParser1 p=new LineParser1('\t');
		ByteFile input=ByteFile.makeByteFile(required(o, "roster"), false);
		try{
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				if(line.length==0 || line[0]=='#'){continue;}
				p.set(line); require(p.terms()==7, "Clean roster must have seven columns");
				final int tid=(int)number(p, 5, Integer.MAX_VALUE);
				final String unit=p.parseString(0), domain=p.parseString(6);
				require(tid>0 && (domain.equals("Bacteria") || domain.equals("Archaea")), "Invalid roster taxonomy");
				require(!byUnit.containsKey(unit) && byTid.get(tid)<0, "Duplicate roster identity");
				final Genome genome=new Genome(tid, domain);
				byUnit.put(unit, genome); byTid.put(tid, genomes.size()); genomes.add(genome);
			}
		}finally{close(input);}
		require(genomes.size()==Integer.parseInt(required(o, "expectedunits")) && !genomes.isEmpty(), "Clean roster size differs");
		input=ByteFile.makeByteFile(required(o, "labels"), false);
		try{
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				if(line.length==0 || line[0]=='#'){continue;}
				p.set(line); require(p.terms()==14, "QuickClade labels require fourteen columns");
				final Genome genome=byUnit.get(p.parseString(0));
				if(genome==null){continue;}
				require(!genome.labeled, "Duplicate eligible QuickClade label"); genome.labeled=true;
				final String status=p.parseString(8);
				require(status.equals("classified") || status.equals("partial") || status.equals("unknown"), "Invalid QuickClade status");
				if(!status.equals("unknown")){require(p.termEquals(genome.domain, 10), "Roster/QuickClade domain mismatch");}
				if(status.equals("classified")){
					genome.phylum=p.parseString(11);
					require(!genome.phylum.isEmpty() && !genome.phylum.equals("unknown"), "Classified label lacks phylum");
				}
			}
		}finally{close(input);}
		final TreeMap<String,Population> groups=new TreeMap<String,Population>();
		for(Genome genome:genomes){
			require(genome.labeled, "Eligible genome has no QuickClade label: "+genome.tid);
			genome.global=population(groups, "global", "-", families);
			genome.domainPop=population(groups, "domain", genome.domain, families);
			genome.phylumPop=genome.phylum==null ? null : population(groups, "phylum", genome.phylum, families);
			genome.global.genomes++; genome.domainPop.genomes++;
			if(genome.phylumPop!=null){genome.phylumPop.genomes++;}
		}
		input=ByteFile.makeByteFile(required(o, "perorg"), false);
		boolean header=false;
		try{
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				if(line.length==0){continue;}
				if(line[0]=='#'){
					require(!header && new String(line, StandardCharsets.US_ASCII).equals(PERORG_HEADER), "Unexpected per-organism header");
					header=true; continue;
				}
				require(header, "Per-organism data precedes header");
				p.set(line); require(p.terms()==16, "Per-organism cache requires sixteen fields");
				final int index=byTid.get((int)number(p, 0, Integer.MAX_VALUE));
				if(index<0){continue;}
				final Genome genome=genomes.get(index);
				require(!genome.counted && p.termEquals(genome.domain, 1), "Duplicate or mismatched per-organism row");
				genome.counted=true;
				p.setBounds(15); final int end=p.b(); int start=p.a(), previous=-1;
				if(end-start==1 && line[start]=='-'){continue;}
				// PerOrgRollup writes an empty final field for genomes with no tracked families;
				// Schema7ExpectedCopyInputsBuilder also retains those genomes in denominators.
				if(start==end){continue;}
				while(start<end){
					int colon=start; while(colon<end && line[colon]!=':' && line[colon]!=';'){colon++;}
					require(colon<end && line[colon]==':', "Malformed family count");
					int stop=colon+1; while(stop<end && line[stop]!=';'){stop++;}
					final int rank=(int)number(line, start, colon, Integer.MAX_VALUE);
					final long copies=number(line, colon+1, stop, Integer.MAX_VALUE);
					require(rank>previous && rank<families && copies>0, "Family ranks must be ordered, unique and positive");
					genome.global.add(rank, copies); genome.domainPop.add(rank, copies);
					if(genome.phylumPop!=null){genome.phylumPop.add(rank, copies);}
					previous=rank; start=stop+1;
				}
			}
		}finally{close(input);}
		for(Genome genome:genomes){require(genome.counted, "Clean genome absent from per-organism cache: "+genome.tid);}
		final String[] inputs=new String[INPUTS.length];
		for(int i=0; i<inputs.length; i++){inputs[i]=required(o, INPUTS[i]);}
		try(OutputFileSetTransaction tx=new OutputFileSetTransaction(new String[]{required(o, "out")}, inputs, false)){
			final ByteStreamWriter writer=new ByteStreamWriter(tx.temporaryPath(0), false, false, false); writer.start();
			try{
				writer.print(new ByteBuilder("#schema\t").append(SCHEMA).append("\n#families\t").append(families)
					.append("\n#family_sha80\t").append(required(o, "familiessha80")).nl());
				for(String name:INPUTS){writer.print(new ByteBuilder("#source\t").append(name).tab().append(required(o, name)).tab().append(required(o, name+"sha80")).nl());}
				for(Population group:groups.values()){
					writer.print(new ByteBuilder("#population\t").append(group.kind).tab().append(group.name).tab().append(group.genomes).nl());
					final ByteBuilder rows=new ByteBuilder();
					for(int rank=0; rank<families; rank++){
						group.validate(rank);
						rows.append(rank).tab().append(group.carriers[rank]).tab().append(group.singletons[rank]).tab().append(group.copies[rank]).nl();
					}
					writer.print(rows);
				}
			}finally{require(!writer.poisonAndWait(), "Requiredness resource write failed");}
			tx.commit();
		}
		System.err.println("REQUIREDNESS_TABLE_PASS genomes="+genomes.size()+" families="+families+" populations="+groups.size());
	}

	/** Loads one frozen resource; population lookup never consults truth labels. */
	static MagQCRequirednessTable load(String path, String pin, int families, String familyPin){
		require(families>0, "Requiredness family width must be positive");
		DigestSuffix.requireSuffix(familyPin, "families");
		require(DigestSuffix.file(path).equals(DigestSuffix.requireSuffix(pin, "requiredness")), "Requiredness resource pin mismatch");
		final TreeMap<String,Population> groups=new TreeMap<String,Population>();
		final ByteFile input=ByteFile.makeByteFile(path, false); final LineParser1 p=new LineParser1('\t');
		Population current=null; int rank=0; boolean schema=false, width=false, family=false;
		try{
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				require(line.length>0, "Blank requiredness resource row");
				p.set(line);
				if(p.termEquals("#schema", 0)){require(!schema && p.terms()==2 && p.termEquals(SCHEMA, 1), "Requiredness schema mismatch"); schema=true;}
				else if(p.termEquals("#families", 0)){require(!width && p.terms()==2 && number(p, 1, Integer.MAX_VALUE)==families, "Requiredness width mismatch"); width=true;}
				else if(p.termEquals("#family_sha80", 0)){require(!family && p.terms()==2 && p.termEquals(familyPin, 1), "Requiredness family order mismatch"); family=true;}
				else if(p.termEquals("#source", 0)){require(p.terms()==4, "Malformed requiredness provenance"); DigestSuffix.requireSuffix(p.parseString(3), "source");}
				else if(p.termEquals("#population", 0)){
					require(schema && width && family && (current==null || rank==families), "Incomplete requiredness population");
					require(p.terms()==4, "Malformed requiredness population header");
					final String kind=p.parseString(1), name=p.parseString(2), key=kind+"\t"+name;
					require(!groups.containsKey(key), "Duplicate requiredness population");
					current=new Population(kind, name, families); current.genomes=number(p, 3, Long.MAX_VALUE);
					require(current.genomes>0, "Empty requiredness population"); groups.put(key, current); rank=0;
				}else{
					require(current!=null && rank<families && p.terms()==4 && number(p, 0, Integer.MAX_VALUE)==rank, "Requiredness ranks differ");
					current.carriers[rank]=number(p, 1, Long.MAX_VALUE); current.singletons[rank]=number(p, 2, Long.MAX_VALUE);
					current.copies[rank]=number(p, 3, Long.MAX_VALUE); current.validate(rank++);
				}
			}
		}finally{close(input);}
		require(current!=null && rank==families, "Incomplete requiredness table");
		return new MagQCRequirednessTable(families, groups);
	}

	/** Keeps the table immutable to consumers; fixtures may construct small exact-count populations. */
	MagQCRequirednessTable(int families_, TreeMap<String,Population> groups_){
		require(families_>0 && groups_!=null, "Requiredness table needs positive width and populations");
		families=families_; groups=new TreeMap<String,Population>(groups_);
		global=requiredPopulation("global", "-"); bacteria=requiredPopulation("domain", "Bacteria"); archaea=requiredPopulation("domain", "Archaea");
		require(global.genomes>=bacteria.genomes && global.genomes-bacteria.genomes==archaea.genomes,
			"Domain denominators must partition the complete corpus");
		for(Population population:groups.values()){
			require(population.carriers.length==families, "Requiredness population width differs");
			for(int rank=0; rank<families; rank++){population.validate(rank);}
		}
		for(int rank=0; rank<families; rank++){
			require(partition(global.carriers[rank], bacteria.carriers[rank], archaea.carriers[rank])
				&& partition(global.singletons[rank], bacteria.singletons[rank], archaea.singletons[rank])
				&& partition(global.copies[rank], bacteria.copies[rank], archaea.copies[rank]),
				"Domain counts do not partition the corpus at rank "+rank);
		}
	}

	/** Unknown/unrepresented predicted phyla use the universal statistic, never a reference TID. */
	Population phylum(String name){final Population p=groups.get("phylum\t"+name); return p==null ? global : p;}
	/** Enumerates the available QuickClade phyla only during feature-resource initialization. */
	java.util.ArrayList<String> phylumNames(){
		final java.util.ArrayList<String> names=new java.util.ArrayList<String>();
		for(Population population:groups.values()){if(population.kind.equals("phylum")){names.add(population.name);}}
		return names;
	}
	/** Retrieves a mandatory domain/global population without synthesizing a denominator. */
	private Population requiredPopulation(String kind, String name){
		final Population p=groups.get(kind+"\t"+name);
		require(p!=null, "Missing requiredness population "+kind+"/"+name); return p;
	}
	/** Checks additive conservation without overflowing while validating untrusted resource counts. */
	private static boolean partition(long total, long first, long second){return total>=first && total-first==second;}

	/** Counts for one frozen population; absent families remain explicit zero rows. */
	static final class Population {
		/** Allocates one complete family-order block for a named population. */
		Population(String kind_, String name_, int families){
			require(families>0 && (kind_.equals("global") || kind_.equals("domain") || kind_.equals("phylum")) && !name_.isEmpty(), "Invalid population identity");
			kind=kind_; name=name_; carriers=new long[families]; singletons=new long[families]; copies=new long[families];
		}
		/** Adds one genome's positive count exactly once after the roster join. */
		void add(int rank, long count){
			assert(count>0 && rank>=0 && rank<carriers.length) : "Sparse cache records contain positive counts at bound family ranks";
			carriers[rank]=Math.addExact(carriers[rank], 1);
			if(count==1){singletons[rank]=Math.addExact(singletons[rank], 1);}
			copies[rank]=Math.addExact(copies[rank], count);
		}
		/** Each non-singleton carrier requires at least two copies; all counts must remain nonnegative. */
		void validate(int rank){
			assert(rank>=0 && rank<carriers.length) : "Population ranks must use the resource's complete family ordering";
			final long c=carriers[rank], s=singletons[rank], n=copies[rank];
			require(genomes>0 && s>=0 && s<=c && c<=genomes && n>=c && n-c>=c-s && (c>0 || n==0),
				"Inconsistent requiredness counts at rank "+rank+" in "+kind+"/"+name);
		}
		/** Presence denominator includes eligible genomes with no member of this family. */
		double presence(int rank){return carriers[rank]/(double)genomes;}
		/** Carrier mean excludes absent genomes, as required by the multicopy feature definition. */
		double carrierMean(int rank){return carriers[rank]==0 ? 0 : copies[rank]/(double)carriers[rank];}
		final String kind, name; long genomes;
		final long[] carriers, singletons, copies;
	}

	/** Shares one mutable accumulator per population during the single corpus pass. */
	private static Population population(TreeMap<String,Population> groups, String kind, String name, int families){
		final String key=kind+"\t"+name; Population p=groups.get(key);
		if(p==null){p=new Population(kind, name, families); groups.put(key, p);} return p;
	}
	/** Validates the existing contiguous zero-based family rank contract. */
	private static int familyCount(String path){
		final ByteFile input=ByteFile.makeByteFile(path, false); final LineParser1 p=new LineParser1('\t'); int count=0;
		try{for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
			if(line.length==0 || line[0]=='#'){continue;} p.set(line);
			require(p.terms()==3 && number(p, 0, Integer.MAX_VALUE)==count, "Family list ranks must be consecutive"); count++;
		}}finally{close(input);} require(count>0, "Empty family list"); return count;
	}
	/** Parses one integer field directly from the reused byte parser. */
	private static long number(LineParser1 p, int field, long max){p.setBounds(field); return number(p.line(), p.a(), p.b(), max);}
	/** Parses unsigned decimal counts, rejecting malformed text and overflow in production too. */
	private static long number(byte[] line, int start, int stop, long max){
		assert(start>=0 && stop<=line.length && max>=0) : "Count ranges must remain within the parser's source row";
		require(start<stop, "Empty count"); long value=0;
		for(int i=start; i<stop; i++){final int d=line[i]-'0'; require(d>=0 && d<10 && value<=(max-d)/10, "Invalid or overflowing count"); value=value*10+d;} return value;
	}
	/** Recognizes the four explicitly bound inputs and their hashes. */
	private static boolean isInput(String key){for(String name:INPUTS){if(key.equals(name) || key.equals(name+"sha80")){return true;}} return false;}
	/** Retrieves a required CLI field without defaulting to an unpinned resource. */
	private static String required(HashMap<String,String> o, String key){final String value=o.get(key); require(value!=null, "Missing "+key); return value;}
	/** Turns delayed reader errors into failures before publishing the output. */
	private static void close(ByteFile input){require(!input.close(), "Requiredness input read failed");}
	/** External resource validation must remain active when assertions are disabled. */
	private static void require(boolean valid, String message){if(!valid){throw new IllegalArgumentException(message);}}
	/** One eligible genome's join state; family observations are streamed rather than retained. */
	private static final class Genome {
		Genome(int tid_, String domain_){tid=tid_; domain=domain_;}
		final int tid; final String domain; String phylum; boolean labeled, counted;
		Population global, domainPop, phylumPop;
	}
	final int families;
	final Population global, bacteria, archaea;
	private final TreeMap<String,Population> groups;
	static final String SCHEMA="magqc_requiredness_v1";
	private static final String[] INPUTS={"perorg", "roster", "labels", "families"};
	private static final String PERORG_HEADER="#tid\tdomain\tlength\tgc\tacgt\tcds\tmapped\tglenSum\tglenSq\tcoding\tr16\tr23\tr5\trother\ttrna\tfamilies";
}
