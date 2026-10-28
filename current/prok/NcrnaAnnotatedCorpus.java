package prok;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicInteger;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import shared.Tools;
import stream.Read;
import structures.ByteBuilder;
import structures.IntList;
import template.ThreadWaiter;

/** Extracts a named Rfam family's complete annotated population using the survey's
 * Parent/interval reconstruction. Source flags are evidence, not independent truth.
 * @author Brian Bushnell, Raiden */
public final class NcrnaAnnotatedCorpus {

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaAnnotatedCorpus.class, false);
		String roster=null, out=null, family=null;int expected=-1, limit=0, threads=8, flank=200, expectedLoci=-1;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0, "Expected flag=value: "+arg);
			final String key=arg.substring(0, eq).toLowerCase(Locale.ROOT), value=arg.substring(eq+1);
			if(key.equals("roster")){roster=value;}else if(key.equals("outprefix")){out=value;}
			else if(key.equals("family")){family=value;}else if(key.equals("expectedgenomes")){expected=Parse.parseIntKMG(value);}
			else if(key.equals("limit")){limit=Parse.parseIntKMG(value);}else if(key.equals("expectedloci")){expectedLoci=Parse.parseIntKMG(value);}
			else if(key.equals("t") || key.equals("threads")){threads=Parse.parseIntKMG(value);}
			else if(key.equals("flank")){flank=Parse.parseIntKMG(value);}else{throw new IllegalArgumentException(key);}
		}
		require(roster!=null && out!=null && family!=null && family.matches("RF[0-9]{5}")
			&& expected>0 && limit>=0 && threads>0 && flank>=0 && flank<=2000 && expectedLoci>=-1,
			"Declare the frozen genome roster, exact Rfam family, sizes, threads and bounded genomic flank");
		Shared.setThreads(threads);Shared.TRIM_READ_DESCRIPTION=false;Read.TO_UPPER_CASE=true;
		ReadWrite.setZipThreads(1);ReadWrite.USE_UNPIGZ=false;
		final NcrnaAnnotatedCorpus tool=new NcrnaAnnotatedCorpus();
		tool.load(roster, expected, limit);tool.extract(family, flank, threads);tool.write(out, family, expectedLoci, flank);
		Shared.closeStream(pp.outstream);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Source Selection       ----------------*/
	/*--------------------------------------------------------------*/

	private void load(final String path, final int expected, final int limit){
		assert(genomes.isEmpty()):"A corpus is bound to one immutable survey roster";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		final ObjectSet<String> names=new ObjectSet<String>(String.class);int count=0;
		try{
			RrnaResourceIO.header(in, "index\tdomain\tphylum\tgenome\tgff\tcompressed_bytes\tmtime_ms\tfeature_rows\tcds_rows\tretained_layers\tcounted_loci\texact_duplicate_components\tlayers_without_noncoding_branch", path);
			for(byte[] row; (row=in.nextLine())!=null;){
				p.set(row);require(p.terms()==13 && p.parseInt(0)==count, "Survey roster order/schema changed");
				final String name=p.parseString(3), gff=p.parseString(4);
				require(names.add(name) && gff.startsWith("/") && gff.endsWith(".gff.gz"), "Unique source names and absolute paired GFF paths required");
				if(limit==0 || count<limit){genomes.add(new Genome(count, name, p.parseString(1), p.parseString(2), gff));}
				count++;
			}
		}finally{require(!in.close(), "Survey roster read failed: "+path);}
		require(count==expected && !genomes.isEmpty(), "Corpus membership must match the complete frozen survey, even in a pilot");
	}

	private void extract(final String family, final int flank, final int threads){
		assert(!genomes.isEmpty()):"Selection must precede genome workers";
		final AtomicInteger next=new AtomicInteger();final ArrayList<Worker> workers=new ArrayList<Worker>();
		for(int i=0; i<Math.min(threads, genomes.size()); i++){workers.add(new Worker(genomes, next, family, flank));}
		require(ThreadWaiter.startAndWait(workers), "All genome workers must finish before corpus publication");
		for(Worker worker:workers){
			if(worker.failure!=null){throw new IllegalStateException("Annotated corpus extraction failed", worker.failure);}
			require(worker.success, "An extraction worker terminated without completing its assigned genomes");
		}
	}

	private static void extractGenome(final Genome g, final String family, final int flank){
		assert(g.loci.isEmpty()):"A source genome must be extracted exactly once";
		final long[] counts=new long[6];
		for(NcrnaFamilySurvey.Locus locus:NcrnaFamilySurvey.readLoci(g.gff, counts)){
			if(locus.family.equals(family)){g.loci.add(new Gene(locus));}
		}
		if(g.loci.isEmpty()){return;}
		require(new File(g.fasta).isFile(), "Missing paired genome: "+g.fasta);
		Collections.sort(g.loci, (a,b)->{
			final int c=a.source.contig.compareTo(b.source.contig);if(c!=0){return c;}
			final int s=a.source.strand.compareTo(b.source.strand);return s!=0 ? s : a.source.spans.compareTo(b.source.spans);
		});
		final ObjectMap<String, Genes> targets=new ObjectMap<String, Genes>(String.class, Genes.class);
		for(Gene gene:g.loci){
			final String key=gene.source.contig;Genes list=targets.get(key);
			if(list==null){list=new Genes();targets.put(key, list);}list.add(gene);
		}
		RrnaResourceIO.stream(g.fasta, r->{
			g.bases+=r.length();g.contigs++;
			final Genes list=targets.get(RrnaResourceIO.first(r.id));if(list==null){return;}
			for(Gene gene:list){
				require(gene.exact==null && gene.stop<=r.length(), "A locus must bind once within its original contig: "+g.name+"/"+gene.source.spans);
				gene.contigLength=r.length();gene.edge=gene.start==1 || gene.stop==r.length();
				gene.exact=new byte[gene.source.length];int offset=0;
				for(int i=0; i<gene.intervals.length; i+=2){
					final int a=gene.intervals[i]-1, length=gene.intervals[i+1]-a;
					System.arraycopy(r.bases, a, gene.exact, offset, length);offset+=length;
				}
				require(offset==gene.exact.length, "Extracted interval union must equal the survey length");
				Tools.toUpperCase(gene.exact);
				final int a=Math.max(0, gene.start-1-flank), z=(int)Math.min(r.length(), (long)gene.stop+flank);
				gene.context=Arrays.copyOfRange(r.bases, a, z);Tools.toUpperCase(gene.context);
				gene.left=gene.start-1-a;gene.right=z-gene.stop;
				if(gene.source.strand.equals("-")){
					AminoAcid.reverseComplementBasesInPlace(gene.exact);AminoAcid.reverseComplementBasesInPlace(gene.context);
					final int left=gene.left;gene.left=gene.right;gene.right=left;
				}
				for(byte base:gene.exact){if(base!='A' && base!='C' && base!='G' && base!='T'){gene.undefined++;}}
				if(gene.intervals.length==2){
					require(gene.context.length-gene.left-gene.right==gene.exact.length,
						"Sense-oriented real flanks must leave the complete annotated interval");
					for(int i=0; i<gene.exact.length; i++){require(gene.exact[i]==gene.context[i+gene.left], "Exact and de-flanked context sequences disagree");}
				}
			}
		});
		for(Gene gene:g.loci){require(gene.exact!=null, "Annotated source seqid was not found in its paired FASTA: "+g.name+"/"+gene.source.contig);}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Output            ----------------*/
	/*--------------------------------------------------------------*/

	private void write(final String prefix, final String family, final int expectedLoci, final int flank){
		assert(!genomes.isEmpty()):"Output denominators come from the declared genome selection";
		RrnaResourceIO.fresh(prefix, ".all.fa.gz", ".full.fa.gz", ".contexts.fa.gz", ".loci.tsv.gz", ".genomes.tsv", ".summary.tsv");
		final ByteStreamWriter all=RrnaResourceIO.writer(prefix+".all.fa.gz"), full=RrnaResourceIO.writer(prefix+".full.fa.gz");
		final ByteStreamWriter contexts=RrnaResourceIO.writer(prefix+".contexts.fa.gz"), rows=RrnaResourceIO.writer(prefix+".loci.tsv.gz");
		final ByteStreamWriter roster=RrnaResourceIO.writer(prefix+".genomes.tsv"), summary=RrnaResourceIO.writer(prefix+".summary.tsv");
		final ByteBuilder b=new ByteBuilder();int raw=0, eligible=0, partial=0, pseudo=0, multipart=0, positive=0, edges=0, ambiguous=0;long territory=0;
		final IntList rawLengths=new IntList(), fullLengths=new IntList();
		try{
			rows.println("locus_id\tgenome\tassembly\ttaxid\tdomain\tphylum\tseqid\tstrand\tintervals1\tlength\tsource_line\tsource_ids\tpartial\tpseudo\tstatus\tcontig_length\tedge\tundefined_bases\tlflank\trflank\tsource_lines");
			roster.println("index\tgenome\tassembly\ttaxid\tdomain\tphylum\tgff\tfasta\tloci\tgenome_bases\tcontigs");
			for(Genome g:genomes){
				if(!g.loci.isEmpty()){positive++;}
				roster.print(b.clear().append(g.index).tab().append(g.name).tab().append(g.assembly).tab().append(g.taxid)
					.tab().append(g.domain).tab().append(g.phylum).tab().append(g.gff).tab().append(g.fasta).tab().append(g.loci.size()).tab().append(g.bases).tab().append(g.contigs).nl());
				int ordinal=0;
				for(Gene gene:g.loci){
					final NcrnaFamilySurvey.Locus l=gene.source;
					final String id=family+"_"+g.index+"_"+ordinal++, status=gene.status();
					raw++;territory+=l.length;rawLengths.add(l.length);if(gene.edge){edges++;}if(gene.undefined>0){ambiguous++;}
					if(status.equals("PSEUDOGENE")){pseudo++;}else if(status.equals("SOURCE_PARTIAL")){partial++;}
					else if(status.equals("MULTIPART_UNRESOLVED")){multipart++;}else{eligible++;}
					b.clear().append('>').append(id).append(" source=").append(g.name).append(" family=").append(family)
						.append(" status=").append(status).nl().append(gene.exact).nl();all.print(b);
					if(status.equals("ANNOTATION_FULL_COPY")){full.print(b);fullLengths.add(l.length);}
					contexts.print(b.clear().append('>').append(id).append(" source=").append(g.name).append(" contig=").append(l.contig)
						.append(" strand=").append(l.strand).append(" start=").append(gene.start).append(" stop=").append(gene.stop)
						.append(" lflank=").append(gene.left).append(" rflank=").append(gene.right).nl().append(gene.context).nl());
					rows.print(b.clear().append(id).tab().append(g.name).tab().append(g.assembly).tab().append(g.taxid).tab().append(g.domain).tab().append(g.phylum)
						.tab().append(l.contig).tab().append(l.strand).tab().append(l.spans).tab().append(l.length).tab().append(l.earliestSourceLine).tab().append(l.sourceIds)
						.tab().append(l.partial).tab().append(l.pseudo).tab().append(status).tab().append(gene.contigLength).tab().append(gene.edge)
						.tab().append(gene.undefined).tab().append(gene.left).tab().append(gene.right).tab().append(l.sourceLines).nl());
				}
			}
			require(raw==eligible+partial+pseudo+multipart, "Every reconstructed family locus needs one explicit disposition");
			require(expectedLoci<0 || raw==expectedLoci, "Extracted raw loci must reconcile with the survey: "+raw+" expected="+expectedLoci);
			summary.println("genomes\tpositive_genomes\tfamily\traw_loci\tterritory_bp\tannotation_full\tsource_partial\tpseudogene\tmultipart_unresolved\tedge_loci\tambiguous_loci\tflank\traw_min\traw_p10\traw_median\traw_p90\traw_max\tfull_min\tfull_p10\tfull_median\tfull_p90\tfull_max");
			b.clear().append(genomes.size()).tab().append(positive).tab().append(family).tab().append(raw).tab().append(territory)
				.tab().append(eligible).tab().append(partial).tab().append(pseudo).tab().append(multipart).tab().append(edges).tab().append(ambiguous).tab().append(flank);
			appendLengths(b, rawLengths);appendLengths(b, fullLengths);summary.print(b.nl());
		}finally{RrnaResourceIO.finish(all, full, contexts, rows, roster, summary);}
		System.err.println("ANNOTATED_CORPUS_PASS genomes="+genomes.size()+" raw="+raw+" annotation_full="+eligible+" territory="+territory);
	}
	/** P10/P90 use floor((n-1)*p); median averages the two middle lengths. */
	private static void appendLengths(final ByteBuilder b, final IntList lengths){
		if(lengths.size==0){b.append("\tNA\tNA\tNA\tNA\tNA");return;}
		lengths.sort();final int n=lengths.size;
		assert(lengths.array[0]>0):"Survey interval unions must have positive length before extraction statistics";
		b.tab().append(lengths.array[0]).tab().append(lengths.array[(n-1)/10])
			.tab().append(((long)lengths.array[(n-1)/2]+lengths.array[n/2])/2.0, 1)
			.tab().append(lengths.array[(int)((n-1)*9L/10)]).tab().append(lengths.array[n-1]);
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Worker extends Thread{
		Worker(final ArrayList<Genome> gs, final AtomicInteger n, final String f, final int pad){genomes=gs;next=n;family=f;flank=pad;}
		@Override
		public void run(){
			try{
				for(int i; (i=next.getAndIncrement())<genomes.size();){
					extractGenome(genomes.get(i), family, flank);
					if((i+1)%1000==0){System.err.println("CORPUS_PROGRESS genomes="+(i+1)+" total="+genomes.size());}
				}
				success=true;
			}catch(Throwable t){failure=t;next.set(genomes.size());}
		}
		final ArrayList<Genome> genomes;final AtomicInteger next;final String family;final int flank;
		Throwable failure;boolean success;
	}
	private static final class Genome{
		Genome(final int i, final String n, final String d, final String p, final String g){
			index=i;name=n;domain=d;phylum=p;gff=g;fasta=g.substring(0, g.length()-7)+".fna.gz";
			final int begin=Math.max(n.lastIndexOf("GCF_"), n.lastIndexOf("GCA_")), version=n.lastIndexOf('_');
			require(begin>=0 && version>begin+4 && n.startsWith("tid_"), "Source identity requires a taxid and assembly accession: "+n);
			assembly=n.substring(begin, version)+"."+n.substring(version+1);
			taxid=Integer.parseInt(n.substring(4, n.indexOf('_', 4)));
		}
		final int index, taxid;final String name, domain, phylum, gff, fasta, assembly;
		final ArrayList<Gene> loci=new ArrayList<Gene>();long bases;int contigs;
	}
	private static final class Gene{
		Gene(final NcrnaFamilySurvey.Locus l){
			source=l;require(l.strand.equals("+") || l.strand.equals("-"), "Sense extraction requires a known source strand");
			final String[] spans=l.spans.split(",");intervals=new int[spans.length*2];
			for(int i=0; i<spans.length; i++){
				final int dash=spans[i].indexOf('-');require(dash>0, "Survey intervals must be positive inclusive ranges");
				intervals[2*i]=Integer.parseInt(spans[i].substring(0, dash));intervals[2*i+1]=Integer.parseInt(spans[i].substring(dash+1));
			}
			start=intervals[0];stop=intervals[intervals.length-1];
		}
		String status(){return source.pseudo ? "PSEUDOGENE" : source.partial ? "SOURCE_PARTIAL" : intervals.length>2 ? "MULTIPART_UNRESOLVED" : "ANNOTATION_FULL_COPY";}
		final NcrnaFamilySurvey.Locus source;final int[] intervals;final int start, stop;
		byte[] exact, context;int contigLength, left, right, undefined;boolean edge;
	}
	private static final class Genes extends ArrayList<Gene>{private static final long serialVersionUID=1L;}
	private static void require(final boolean ok, final String why){if(!ok){throw new IllegalArgumentException(why);}}
	private final ArrayList<Genome> genomes=new ArrayList<Genome>();
}
