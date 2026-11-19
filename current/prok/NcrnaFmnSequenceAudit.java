package prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import map.LongIntMap;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import stream.Read;
import structures.ByteBuilder;
import structures.IntList;

/** Exact saved-locus extraction and both-strand31mer presence audit.
 * Native sequence owners read original FASTAs, never a substituted CDS/window file.
 * A shared31mer is a seed, not a biological annotation. The reported best length
 * is one contiguous exact normalized DNA match to one query, extended from a seed.
 * No alignment, CM search, split change or resource promotion is performed.
 * @author Ganyu
 */
public final class NcrnaFmnSequenceAudit{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaFmnSequenceAudit.class, false);
		try{
			final Config c=new Config();
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String a=pair[0].toLowerCase(Locale.ROOT);
				final String b=pair.length>1 ? pair[1] : null;
				if(a.equals("cores")){c.cores=b;}
				else if(a.equals("metadata")){c.metadata=b;}
				else if(a.equals("census")){c.census=b;}
				else if(a.equals("previous")){c.previous=b;}
				else if(a.equals("expectedqueries")){c.queries=Parse.parseIntKMG(b);}
				else if(a.equals("expecteddiscovery")){c.discovery=Parse.parseIntKMG(b);}
				else{require(parser.parse(arg, a, b), "Unknown FMN sequence audit option: "+arg);}
			}
			c.lineages=parser.in1;
			c.output=parser.out1;
			require(parser.in2==null && parser.out2==null && !parser.overwrite && !parser.append
				&& parser.maxReads<0 && parser.samplerate==1, "Audit complete unpaired populations without sampling or overwrite");
			run(c);
			pp.outstream.println("NCRNA_FMN_SEQUENCE_AUDIT_COMPLETE_NOT_CM_OR_SPLIT");
		}finally{Shared.closeStream(pp.outstream);}
	}

	public static void run(final Config c){
		require(c!=null && c.lineages!=null && c.cores!=null && c.metadata!=null && c.census!=null
			&& c.previous!=null && c.output!=null && c.queries>0 && c.discovery>0, "Declare the complete saved inputs and positive expected populations");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Sequence audit requires native assertions");
		RrnaResourceIO.fresh(c.output, ".selection.tsv", ".queries.fa", ".presence.tsv", ".population.tsv",
			".metadata.tsv.gz", ".discovery.census.tsv", ".expected.tsv", ".summary.tsv");
		final NcrnaFmnSequenceAudit audit=new NcrnaFmnSequenceAudit();
		audit.loadLineages(c.lineages);
		audit.metadata(c.metadata, c.output);
		audit.select(c.queries, c.discovery);
		audit.census(c.census, c.output);
		final int prior=audit.previous(c.previous);
		audit.cores(c.cores, c.output, c.queries);
		final boolean previousUpper=Read.TO_UPPER_CASE;
		final boolean previousMask=Read.LOWER_CASE_TO_N;
		try{
			Read.TO_UPPER_CASE=false;
			Read.LOWER_CASE_TO_N=false;
			audit.extract(c.output);
			audit.match(c.output, prior);
		}finally{
			Read.TO_UPPER_CASE=previousUpper;
			Read.LOWER_CASE_TO_N=previousMask;
		}
	}

	private void loadLineages(final String path){
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, LINEAGE_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 11);
				final Genome g=new Genome();
				g.cohort=text(0);
				g.id=p.parseInt(1);
				g.path=text(2);
				g.group=p.parseInt(3);
				g.assembly=text(4);
				g.genus=text(8);
				g.cores=p.parseInt(10);
				require((g.cohort.equals("evaluation") || g.cohort.equals("discovery")) && g.id>0 && g.group>=0
					&& g.path.startsWith("/") && g.cores>=0 && p.termEquals("GTDB_ACCESSION_MATCH", 9),
					"Genus selection requires complete original accession labels and intact input identities");
				require(byId.get(g.key())==null, "Duplicate lineage input: "+g.key());
				byId.put(g.key(), g);
				all.add(g);
			}
		}
	}

	private void metadata(final String path, final String out){
		final ObjectMap<NcrnaReferenceGeometry.Target, String[]> wanted=
			new ObjectMap<NcrnaReferenceGeometry.Target, String[]>(NcrnaReferenceGeometry.Target.class, String[].class);
		final NcrnaReferenceGeometry.Target probe=new NcrnaReferenceGeometry.Target();
		for(final Genome g:all){
			final byte[] accession=utf8(g.assembly);
			probe.set(accession, 0, accession.length);
			wanted.put(probe.copy(), null);
		}
		try(Input in=new Input(path)){
			final byte[] header=in.reader.nextLine();
			require(header!=null, "GTDB metadata needs its original header");
			p.set(header);
			final int terms=p.terms();
			final int[] columns=new int[META.length];
			Arrays.fill(columns, -1);
			for(int i=0; i<terms; i++){
				for(int j=0; j<META.length; j++){
					if(p.termEquals(META[j], i)){
						require(columns[j]<0, "Repeated metadata column");
						columns[j]=i;
					}
				}
			}
			for(final int col:columns){require(col>=0, "Required original GTDB metadata column missing");}
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, terms);
				p.setBounds(columns[0]);
				int start=p.a();
				if(p.b()-start>3 && row[start+2]=='_' && ((row[start]=='R' && row[start+1]=='S') || (row[start]=='G' && row[start+1]=='B'))){start+=3;}
				probe.set(row, start, p.b());
				if(!wanted.contains(probe)){continue;}
				require(wanted.get(probe)==null, "Duplicate original metadata accession");
				final String[] values=new String[META.length];
				for(int j=0; j<columns.length; j++){values[j]=text(columns[j]);}
				wanted.put(probe.copy(), values);
			}
		}
		final ObjectMap<String, long[]> bins=new ObjectMap<String, long[]>(String.class, long[].class);
		final ArrayList<String> keys=new ArrayList<String>();
		try(Output detail=new Output(out+".metadata.tsv.gz")){
			b.clear().append("cohort\tinput_id\tassembly\tgenus");
			for(final String field:META){b.tab().append(field);}
			detail.writer.print(b.nl());
			for(final Genome g:all){
				final byte[] accession=utf8(g.assembly);
				probe.set(accession, 0, accession.length);
				g.meta=wanted.get(probe);
				require(g.meta!=null, "Original GTDB metadata missing accession: "+g.assembly);
				b.clear().append(g.cohort).tab().append(g.id).tab().append(g.assembly).tab().append(utf8(g.genus));
				for(final String value:g.meta){b.tab().append(utf8(value));}
				detail.writer.print(b.nl());
				for(int field=1; field<=4; field++){
					final String key=g.cohort+'\t'+META[field]+'\t'+g.meta[field];
					long[] count=bins.get(key);
					if(count==null){
						count=new long[1];
						bins.put(key, count);
						keys.add(key);
					}
					count[0]++;
				}
			}
		}
		java.util.Collections.sort(keys);
		try(Output table=new Output(out+".population.tsv")){
			table.writer.println("cohort\tmetadata_field\toriginal_value\tinput_records");
			for(final String key:keys){table.writer.print(b.clear().append(utf8(key)).tab().append(bins.get(key)[0]).nl());}
		}
	}

	private void select(final int expectedQueries, final int expectedDiscovery){
		final ObjectSet<String> groupKeys=new ObjectSet<String>(String.class);
		int chryseo=0;
		int elizabeth=0;
		int empedo=0;
		for(final Genome g:all){
			final boolean isEvaluation=g.cohort.equals("evaluation");
			if(isEvaluation ? !g.genus.equals("Chryseobacterium") || g.cores==0
				: !g.genus.equals("Chryseobacterium") && !g.genus.equals("Elizabethkingia") && !g.genus.equals("Empedobacter")){continue;}
			require(groupKeys.add(g.cohort+'\t'+g.group), "Selected diagnostic inputs must be distinct assembly groups within each cohort");
			g.selected=true;
			selected.add(g);
			if(isEvaluation){evaluation.add(g);}else{
				discovery.add(g);
				if(g.genus.equals("Chryseobacterium")){chryseo++;}
				else if(g.genus.equals("Elizabethkingia")){elizabeth++;}
				else{empedo++;}
			}
		}
		require(evaluation.size()==expectedQueries && chryseo==expectedDiscovery && elizabeth==1 && empedo==1,
			"The declared positive evaluation and full discovery genus populations changed");
		java.util.Collections.sort(selected, (a, z)->a.key().compareTo(z.key()));
		java.util.Collections.sort(discovery, (a, z)->Integer.compare(a.id, z.id));
	}

	private void census(final String path, final String out){
		try(Input in=new Input(path); Output subset=new Output(out+".discovery.census.tsv")){
			RrnaResourceIO.header(in.reader, NcrnaEvidenceCatalog.CENSUS_HEADER, path);
			subset.writer.println(NcrnaEvidenceCatalog.CENSUS_HEADER);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 5);
				final Genome g=byId.get("discovery\t"+p.parseInt(0));
				if(g==null || !g.selected){continue;}
				require(g.expectedBases<0 && text(1).equals(g.path), "Original scan census must bind the exact selected file once");
				g.fileBytes=p.parseLong(2);
				g.expectedBases=p.parseLong(3);
				g.expectedContigs=p.parseLong(4);
				require(g.fileBytes>0 && g.expectedBases>0 && g.expectedContigs>0, "Scan dimensions must be positive");
				subset.writer.print(b.clear().append(row).nl());
			}
		}
		for(final Genome g:discovery){require(g.expectedBases>0, "Missing discovery scan census input: "+g.id);}
		try(Output table=new Output(out+".selection.tsv")){
			table.writer.println("cohort\tinput_id\tfasta\tassembly_group\tassembly\tgenus\tCM_admitted_cores");
			for(final Genome g:selected){table.writer.print(b.clear().append(g.cohort).tab().append(g.id).tab().append(utf8(g.path))
				.tab().append(g.group).tab().append(g.assembly).tab().append(utf8(g.genus)).tab().append(g.cores).nl());}
		}
	}

	private int previous(final String path){
		int matches=0;
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, NcrnaConstructionMembers.CROSS_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 12);
				final Genome g=byId.get(text(0)+'\t'+p.parseInt(1));
				require(g!=null && text(2).equals(g.path) && p.parseInt(3)==g.group && text(4).equals(g.assembly),
					"Prior crossed selection must bind the same original inputs and assemblies");
				if(g.selected && g.cohort.equals("discovery")){matches++;}
			}
		}
		return matches;
	}

	private void cores(final String path, final String out, final int expected){
		final ObjectSet<String> ids=new ObjectSet<String>(String.class);
		try(Input in=new Input(path); Output selectedCores=new Output(out+".expected.tsv")){
			RrnaResourceIO.header(in.reader, NcrnaConstructionMembers.EXTRACT_HEADER, path);
			selectedCores.writer.println(NcrnaConstructionMembers.EXTRACT_HEADER);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 21);
				final Genome g=byId.get(text(1)+'\t'+p.parseInt(2));
				if(g==null || !g.selected || !g.cohort.equals("evaluation")){continue;}
				final String admission=text(12);
				if(!token(admission, "CM_FULL_CLEAN")){continue;}
				require(text(3).equals(g.path) && p.parseInt(4)==g.group && p.termEquals("true", 14), "Selected core must bind the native-intact source assembly");
				final Query q=new Query();
				q.id=text(0);
				q.seqid=text(8);
				q.minus=p.termEquals("-", 9);
				require(q.minus || p.termEquals("+", 9), "Core strand must be explicit");
				final String interval=text(10);
				final int dash=interval.indexOf('-');
				require(dash>0 && interval.indexOf(',')<0, "Native full FMN probes require one intact interval");
				q.start=Integer.parseInt(interval.substring(0, dash))-1;
				q.stop=Integer.parseInt(interval.substring(dash+1));
				require(q.start>=0 && q.stop-q.start==p.parseInt(11) && q.stop-q.start>=31 && ids.add(q.id), "Core coordinates/identity must be intact and unique");
				g.queries.add(q);
				queries.add(q);
				selectedCores.writer.print(b.clear().append(row).nl());
			}
		}
		require(queries.size()==expected, "Expected native locus population differs from the declared count");
		for(final Genome g:evaluation){require(g.queries.size()==g.cores, "Selected positive input lost native core metadata: "+g.id);}
	}

	private void extract(final String out){
		for(final Genome g:evaluation){
			final ObjectSet<String> seen=new ObjectSet<String>(String.class);
			NcrnaDdlSavedData.sequences(g.path, read->{
				final String seqid=RrnaResourceIO.first(read.id);
				require(seen.add(seqid), "Duplicate original contig token in extraction input: "+g.path);
				for(final Query q:g.queries){
					if(!q.seqid.equals(seqid)){continue;}
					require(q.bases==null && q.stop<=read.bases.length, "Native locus coordinates exceed original contig");
					q.bases=Arrays.copyOfRange(read.bases, q.start, q.stop);
					for(int i=0; i<q.bases.length; i++){
						final int code=RrnaResourceIO.base(q.bases[i]);
						require(code>=0, "Undefined bases in an exact FMN probe: "+q.id);
						q.bases[i]=DNA[code];
					}
					if(q.minus){q.bases=reverse(q.bases);}
				}
			});
		}
		try(Output fasta=new Output(out+".queries.fa")){
			for(final Query q:queries){
				require(q.bases!=null, "Native locus contig not found: "+q.id);
				fasta.writer.print(b.clear().append('>').append(utf8(q.id)).nl().append(q.bases).nl());
			}
		}
	}

	private void match(final String out, final int previous){
		final byte[][] probes=new byte[queries.size()][];
		for(int i=0; i<probes.length; i++){probes[i]=queries.get(i).bases;}
		final Matcher matcher=new Matcher(probes);
		int hitGenomes=0;
		try(Output result=new Output(out+".presence.tsv")){
			result.writer.println("input_id\tassembly\tgenus\tfasta\tscanned_bases\tscanned_contigs\tundefined_bases\tlowercase_bases\tGC_bases\tGTDB_genome_size\tGTDB_contig_count\tNCBI_total_length\tNCBI_contig_count\tbest_exact_length\tquery_id\tquery_length\tcontig\tstart1\tstop1\tquery_orientation\tinterpretation");
			for(final Genome g:discovery){
				require(new File(g.path).length()==g.fileBytes, "Original scan FASTA file size changed: "+g.path);
				final long[] stats=new long[5];
				final Match best=new Match();
				final ObjectSet<String> seen=new ObjectSet<String>(String.class);
				NcrnaDdlSavedData.sequences(g.path, read->{
					final String token=RrnaResourceIO.first(read.id);
					require(seen.add(token), "Duplicate discovery contig token");
					stats[0]+=read.bases.length;
					stats[1]++;
					for(final byte base:read.bases){
						final int code=RrnaResourceIO.base(base);
						if(code<0){stats[2]++;}
						if(base>='a' && base<='z'){stats[3]++;}
						if(code==1 || code==2){stats[4]++;}
					}
					matcher.scan(read.bases, token, best);
				});
				require(stats[0]==g.expectedBases && stats[1]==g.expectedContigs, "Original scan dimensions changed for "+g.assembly);
				if(best.length>=31){hitGenomes++;}
				b.clear().append(g.id).tab().append(g.assembly).tab().append(utf8(g.genus)).tab().append(utf8(g.path));
				for(final long value:stats){b.tab().append(value);}
				for(int i=5; i<=8; i++){b.tab().append(utf8(g.meta[i]));}
				b.tab().append(best.length).tab().append(best.query<0 ? "NA" : queries.get(best.query).id).tab()
					.append(best.query<0 ? 0 : probes[best.query].length).tab().append(utf8(best.contig)).tab()
					.append(best.start+1).tab().append(best.start+best.length).tab().append(best.minus ? "-" : "+")
					.append("\tEXACT_SEEDED_SUBSTRING_NOT_FMN_ANNOTATION\n");
				result.writer.print(b);
			}
		}
		try(Output summary=new Output(out+".summary.tsv")){
			summary.writer.println("evaluation_queries\tdiscovery_inputs\tdiscovery_with_exact31mer\tprior_crossed_discovery_overlap\tstatus");
			summary.writer.print(b.clear().append(queries.size()).tab().append(discovery.size()).tab().append(hitGenomes).tab().append(previous)
				.append("\tEXACT_SEQUENCE_AUDIT_NOT_CM_OR_SPLIT\n"));
		}
	}

	/** Primitive exact31mer lookup with query-position chains and exact extension.
	 * Opposite query orientations are indexed separately; contigs never join. */
	public static final class Matcher{
		public Matcher(final byte[][] probes){
			require(probes!=null && probes.length>0, "Exact matching requires nonempty probes");
			oriented=new byte[probes.length*2][];
			for(int q=0; q<probes.length; q++){
				require(probes[q]!=null && probes[q].length>=31, "Probe is shorter than the exact31mer seed");
				oriented[2*q]=probes[q].clone();
				oriented[2*q+1]=reverse(probes[q]);
				for(int strand=0; strand<2; strand++){
					final byte[] bases=oriented[2*q+strand];
					long key=0;
					for(int i=0; i<bases.length; i++){
						final int code=RrnaResourceIO.base(bases[i]);
						require(code>=0, "Exact query seeds cannot contain undefined bases");
						key=((key<<2)|code)&MASK;
						if(i<30){continue;}
						final int entry=positions.size;
						next.add(heads.get(key));
						positions.add(i-30);
						owners.add(2*q+strand);
						heads.put(key, entry);
					}
				}
			}
		}
		public void scan(final byte[] bases, final String contig, final Match best){
			require(bases!=null && contig!=null && best!=null, "Each original contig needs separate exact-match bounds");
			long key=0;
			int valid=0;
			for(int i=0; i<bases.length; i++){
				final int code=RrnaResourceIO.base(bases[i]);
				if(code<0){
					key=0;
					valid=0;
					continue;
				}
				key=((key<<2)|code)&MASK;
				if(++valid<31){continue;}
				final int start=i-30;
				for(int entry=heads.get(key); entry>=0; entry=next.get(entry)){
					final int owner=owners.get(entry);
					final byte[] query=oriented[owner];
					final int pos=positions.get(entry);
					// Interior seeds of this same exact run were already extended at its left edge.
					if(start>0 && pos>0 && RrnaResourceIO.base(bases[start-1])==RrnaResourceIO.base(query[pos-1])){continue;}
					int length=31;
					while(start+length<bases.length && pos+length<query.length
						&& RrnaResourceIO.base(bases[start+length])==RrnaResourceIO.base(query[pos+length])){length++;}
					if(length>best.length){
						best.length=length;
						best.query=owner/2;
						best.minus=(owner&1)!=0;
						best.start=start;
						best.contig=contig;
					}
				}
			}
		}
		private final byte[][] oriented;
		private final LongIntMap heads=new LongIntMap();
		private final IntList positions=new IntList();
		private final IntList owners=new IntList();
		private final IntList next=new IntList();
	}
	public static final class Match{
		public int length;
		public int query=-1;
		public int start=-1;
		public boolean minus;
		public String contig="NA";
	}
	public static byte[] reverse(final byte[] bases){
		require(bases!=null, "Reverse complement needs actual sequence bytes");
		final byte[] out=new byte[bases.length];
		for(int i=0; i<bases.length; i++){
			final int code=RrnaResourceIO.base(bases[bases.length-1-i]);
			require(code>=0, "Undefined query bases cannot define an exact reverse complement");
			out[i]=DNA[3-code];
		}
		return out;
	}
	private void parse(final byte[] row, final int fields){
		p.set(row);
		require(p.terms()==fields, "Saved metadata schema changed");
	}
	private String text(final int field){return new String(p.parseByteArray(field), StandardCharsets.UTF_8);}
	private static byte[] utf8(final String text){return text.getBytes(StandardCharsets.UTF_8);}
	private static boolean token(final String values, final String target){
		assert(values!=null && target!=null && !target.isEmpty()) : "Admission matching requires an exact nonempty delimiter-bounded token";
		int start=0;
		while(start<values.length()){
			int stop=values.indexOf(';', start);
			if(stop<0){stop=values.length();}
			if(stop-start==target.length() && values.regionMatches(start, target, 0, target.length())){return true;}
			start=stop+1;
		}
		return false;
	}
	private static void require(final boolean ok, final String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final class Genome{
		String cohort, path, assembly, genus;
		int id, group, cores;
		boolean selected;
		String[] meta;
		long fileBytes, expectedBases=-1, expectedContigs;
		final ArrayList<Query> queries=new ArrayList<Query>();
		String key(){return cohort+'\t'+id;}
	}
	private static final class Query{
		String id, seqid;
		int start, stop;
		boolean minus;
		byte[] bases;
	}
	public static final class Config{
		public String lineages, cores, metadata, census, previous, output;
		public int queries=152, discovery=47;
	}
	private final LineParser1 p=new LineParser1('\t');
	private final ByteBuilder b=new ByteBuilder();
	private final ArrayList<Genome> all=new ArrayList<Genome>();
	private final ArrayList<Genome> selected=new ArrayList<Genome>();
	private final ArrayList<Genome> evaluation=new ArrayList<Genome>();
	private final ArrayList<Genome> discovery=new ArrayList<Genome>();
	private final ArrayList<Query> queries=new ArrayList<Query>();
	private final ObjectMap<String, Genome> byId=new ObjectMap<String, Genome>(String.class, Genome.class);
	private static final byte[] DNA={'A', 'C', 'G', 'T'};
	private static final long MASK=0x3fffffffffffffffL;
	private static final String LINEAGE_HEADER="cohort\tinput_id\tfasta\tassembly_group\tassembly\tclass\torder\tfamily\tgenus\ttaxonomy_status\tCM_admitted_cores";
	private static final String[] META={"accession", "ncbi_genome_category", "ncbi_submitter", "ncbi_assembly_level", "checkm2_completeness",
		"genome_size", "contig_count", "ncbi_total_length", "ncbi_contig_count", "ncbi_genome_representation", "gc_percentage"};
}
