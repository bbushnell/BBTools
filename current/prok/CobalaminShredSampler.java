package prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Locale;
import java.util.Random;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import structures.ByteBuilder;
import structures.IntList;

/** Project-only genomic background sampler; source masking is not a negative label.
 * Balanced genome allocation, empirical lengths, uniform eligible starts within each
 * assigned genome, with replacement. The entire output requires both-strand CM screening.
 * @author Raiden
 */
public final class CobalaminShredSampler {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args, CobalaminShredSampler.class, false);
		final CobalaminShredSampler x=new CobalaminShredSampler();
		String genomes=null, panel=null, lengths=null, loci=null, out=null;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0, "Use flag=value");
			final String k=arg.substring(0, eq).toLowerCase(Locale.ROOT), v=arg.substring(eq+1);
			if(k.equals("genomes")){genomes=v;}else if(k.equals("panel")){panel=v;}else if(k.equals("lengths")){lengths=v;}
			else if(k.equals("loci")){loci=v;}else if(k.equals("outprefix")){out=v;}else if(k.equals("arm")){x.arm=v;}
			else if(k.equals("count")){x.count=Integer.parseInt(v);}else if(k.equals("flank")){x.flank=Integer.parseInt(v);}
			else if(k.equals("genomecount")){x.genomeCount=Integer.parseInt(v);}else if(k.equals("expectedgenomes")){x.expectedGenomes=Integer.parseInt(v);}
			else if(k.equals("expectedpanel")){x.expectedPanel=Integer.parseInt(v);}else if(k.equals("expectedlengths")){x.expectedLengths=Integer.parseInt(v);}
			else if(k.equals("expectedloci")){x.expectedLoci=Integer.parseInt(v);}else if(k.equals("seed")){x.seed=Long.parseLong(v);}
			else{throw new IllegalArgumentException("Unknown option: "+k);}
		}
		require(genomes!=null && panel!=null && lengths!=null && loci!=null && out!=null && x.arm!=null,
			"Bind whole genome roster, panel FASTA sha80 ledger, empirical-length FASTA, three source locus ledgers, output and arm");
		require((x.arm.equals("panel") || x.arm.equals("nonpanel")) && x.count>0 && x.flank>=0 && x.genomeCount>0 && x.count>=x.genomeCount
			&& x.expectedGenomes>=x.genomeCount && x.expectedPanel>0 && x.expectedLengths>0 && x.expectedLoci>0 && x.seed!=Long.MIN_VALUE,
			"Declare positive population/sample counts, nonnegative flank, seed and panel/nonpanel arm");
		Shared.setThreads(1);Shared.TRIM_READ_DESCRIPTION=false;Read.TO_UPPER_CASE=true;
		try{x.roster(genomes, panel);x.lengths(lengths);x.masks(loci);x.write(out);}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	private void roster(String path, String panelPath){
		final ObjectSet<String> panelFiles=new ObjectSet<String>(String.class), panelAssemblies=new ObjectSet<String>(String.class), panelTaxids=new ObjectSet<String>(String.class);
		final LineParser1 p=new LineParser1('\t');ByteFile in=RrnaResourceIO.open(panelPath);
		try{for(byte[] row; (row=in.nextLine())!=null; ){
			p.set(row);require(p.terms()==2 && p.parseString(0).matches("[0-9a-f]{20}") && p.parseString(1).startsWith("/") && panelFiles.add(p.parseString(1)),
				"Panel paths must be unique absolute FASTA identities with sha80 provenance");
		}}finally{require(!in.close(), "Panel ledger read failed");}
		require(panelFiles.size()==expectedPanel, "Panel population changed");
		final ArrayList<Genome> all=new ArrayList<Genome>();final ObjectSet<String> seen=new ObjectSet<String>(String.class), paths=new ObjectSet<String>(String.class);
		in=RrnaResourceIO.open(path);int matched=0;
		try{RrnaResourceIO.header(in, NcrnaCorpusSample.GENOME_HEADER, path);for(byte[] row; (row=in.nextLine())!=null; ){
			p.set(row);require(p.terms()==11, "Expected complete source genome roster schema");
			final Genome g=new Genome(p.parseString(1), p.parseString(2), p.parseString(3), p.parseString(7));
			require(g.fasta.startsWith("/") && seen.add(g.name) && paths.add(g.fasta), "Unique genome names and absolute FASTA paths are required");
			g.panel=panelFiles.contains(g.fasta);all.add(g);
			if(g.panel){matched++;panelAssemblies.add(stem(g.assembly));panelTaxids.add(g.taxid);}
		}}finally{require(!in.close(), "Genome roster read failed");}
		require(all.size()==expectedGenomes && matched==expectedPanel, "Full roster and panel membership must match declared populations");
		for(Genome g:all){
			if(arm.equals("panel")?g.panel:!g.panel && !panelAssemblies.contains(stem(g.assembly)) && !panelTaxids.contains(g.taxid)){selected.add(g);}
		}
		require(selected.size()>=genomeCount && (!arm.equals("panel") || genomeCount==expectedPanel), "Panel arm uses every panel genome; nonpanel arm needs enough source-disjoint genomes");
		eligibleGenomes=selected.size();Collections.sort(selected, (a, b)->a.name.compareTo(b.name));Collections.shuffle(selected, new Random(seed^0x6732a74bL));
		while(selected.size()>genomeCount){selected.remove(selected.size()-1);}
		for(Genome g:selected){byGenome.put(g.name, g);}
	}
	private void lengths(String path){
		final ObjectSet<String> ids=new ObjectSet<String>(String.class);
		RrnaResourceIO.stream(path, r->{
			require(r.bases!=null && r.length()>0 && ids.add(RrnaResourceIO.first(r.id)), "Empirical lengths require unique intact positive records");
			lengths.add(r.length());
		});
		require(lengths.size==expectedLengths, "Empirical positive population changed: "+lengths.size+" expected="+expectedLengths);
		// Sorting makes schedules independent of FASTA row order without changing frequencies.
		Arrays.sort(lengths.array, 0, lengths.size);final Random random=new Random(seed);schedule=new int[count];
		for(int i=0; i<count; i++){schedule[i]=lengths.get(random.nextInt(lengths.size));}
	}
	private void masks(String files){
		final String[] paths=files.split(",");require(paths.length==3, "Supply RF00174, RF01482 and RF01689 source ledgers in that order");
		final LineParser1 p=new LineParser1('\t');final ObjectSet<String> ids=new ObjectSet<String>(String.class);int observed=0;
		for(int f=0; f<3; f++){
			final ByteFile in=RrnaResourceIO.open(paths[f]);
			try{RrnaResourceIO.header(in, NcrnaCmLocusAudit.HEADER, paths[f]);for(byte[] row; (row=in.nextLine())!=null; ){
				p.set(row);require(p.terms()==21, "All source masks need the original annotation schema");
				final String id=p.parseString(0);require(id.startsWith(FAMILIES[f]+"_") && ids.add(id), "Source mask IDs must retain unique declared family provenance");observed++;
				final Genome g=byGenome.get(p.parseString(1));if(g==null){continue;}
				require(p.termEquals(g.assembly, 2) && p.termEquals(g.taxid, 3), "Mask source identity must match selected genome: "+id);
				final String contig=p.parseString(6);Mask mask=g.masks.get(contig);
				if(mask==null){mask=new Mask(p.parseInt(15));g.masks.put(contig, mask);g.maskNames.add(contig);}
				require(mask.length==p.parseInt(15) && mask.length>0, "Source masks must agree on contig length: "+id);
				// Mask the full annotation span, including any internal interval gaps.
				int lo=Integer.MAX_VALUE, hi=-1;
				for(String interval:p.parseString(8).split(",")){
					final int dash=interval.indexOf('-');require(dash>0, "Source interval must contain inclusive bounds: "+id);
					final int a=Integer.parseInt(interval.substring(0, dash)), z=Integer.parseInt(interval.substring(dash+1));
					require(a>0 && z>=a && z<=mask.length, "Source mask must lie within its contig: "+id);lo=Math.min(lo, a);hi=Math.max(hi, z);
				}
				mask.spans.add(new int[]{lo, hi});g.sourceLoci++;
			}}finally{require(!in.close(), "Source mask read failed");}
		}
		require(observed==expectedLoci, "The complete three-family annotation mask population changed");
	}
	private void write(String prefix){
		RrnaResourceIO.fresh(prefix, ".cores.fa.gz", ".contexts.fa.gz", ".shreds.tsv", ".genomes.tsv", ".lengths.tsv", ".summary.tsv");
		final ByteStreamWriter rows=RrnaResourceIO.writer(prefix+".shreds.tsv"), genomes=RrnaResourceIO.writer(prefix+".genomes.tsv"), hist=RrnaResourceIO.writer(prefix+".lengths.tsv"), summary=RrnaResourceIO.writer(prefix+".summary.tsv");
		final Euk18sShredSet.Output core=new Euk18sShredSet.Output(prefix+".cores.fa.gz", 0), context=new Euk18sShredSet.Output(prefix+".contexts.fa.gz", flank);
		final ByteBuilder b=new ByteBuilder();final Random positions=new Random(seed^0x1b93f045L);long coreBases=0, contextBases=0;int written=0;
		try{
			rows.println(HEADER);genomes.println("arm\tgenome\tassembly\ttaxid\tfasta\tsource_mask_loci\tshreds\tbases");
			for(int gi=0; gi<selected.size(); gi++){
				final Genome g=selected.get(gi);final ArrayList<Contig> contigs=new ArrayList<Contig>();final ObjectSet<String> names=new ObjectSet<String>(String.class);
				RrnaResourceIO.stream(g.fasta, r->{
					final String name=RrnaResourceIO.first(r.id);require(r.bases!=null && names.add(name), "Source FASTA needs unique named contigs");
					final Mask source=g.masks.get(name);require(source==null || source.length==r.length(), "Mask coordinates must match actual FASTA contig length: "+name);
					final Euk18sShredSet.Contig mask=new Euk18sShredSet.Contig(name, r.length(), 0, 1, 1);
					if(source!=null){mask.mask.addAll(source.spans);}
					for(int i=0; i<r.length(); ){
						if(RrnaResourceIO.base(r.bases[i])>=0){i++;continue;}
						final int start=i++;while(i<r.length() && RrnaResourceIO.base(r.bases[i])<0){i++;}mask.mask.add(new int[]{start+1, i});
					}
					mask.prepareMask(0);contigs.add(new Contig(name, r.bases, mask));
				});
				for(String name:g.maskNames){require(names.contains(name), "Every selected source annotation contig must occur in its FASTA: "+name);}
				int n=0;long bases=0;
				for(int i=gi; i<count; i+=genomeCount){
					final int length=schedule[i];final ArrayList<Start> starts=startSpace(contigs, length);
					require(!starts.isEmpty(), "Assigned length has no unmasked start in genome; do not redraw its length: "+g.name+" length="+length);
					final long ticket=Euk18sShredSet.bounded(positions, starts.get(starts.size()-1).end);int a=0, z=starts.size();
					while(a<z){final int mid=(a+z)>>>1;if(starts.get(mid).end<=ticket){a=mid+1;}else{z=mid;}}
					final Start s=starts.get(a);final Contig c=s.contig;final int lo=s.start+(int)(ticket-(a==0?0:starts.get(a-1).end)), hi=lo+length-1;
					require(!c.mask.overlaps(lo, hi), "Complement-start sampler must exclude all source annotations and ambiguous core bases");
					final int left=(int)Math.max(1L, (long)lo-flank), right=(int)Math.min(c.bases.length, (long)hi+flank);final String id=arm+"_"+i;
					core.add(new Read(Arrays.copyOfRange(c.bases, lo-1, hi), null, id, i));context.add(new Read(Arrays.copyOfRange(c.bases, left-1, right), null, id, i));
					rows.print(b.clear().append(id).tab().append(arm).tab().append(g.name).tab().append(g.assembly).tab().append(g.taxid).tab().append(c.name)
						.tab().append(lo).tab().append(hi).tab().append(left).tab().append(right).tab().append(lo-left+1).tab().append(hi-left+1).tab().append(length).nl());
					n++;written++;bases+=length;coreBases+=length;contextBases+=right-left+1;
				}
				genomes.print(b.clear().append(arm).tab().append(g.name).tab().append(g.assembly).tab().append(g.taxid).tab().append(g.fasta).tab().append(g.sourceLoci).tab().append(n).tab().append(bases).nl());
			}
			require(written==count, "Every fixed length schedule entry must become one shred");
			final IntList empirical=new IntList(), sampled=new IntList();for(int i=0; i<lengths.size; i++){empirical.increment(lengths.get(i), 1);}for(int n:schedule){sampled.increment(n, 1);}
			hist.println("length\tpositive_corpus\tbackground_sample");for(int i=0; i<Math.max(empirical.size, sampled.size); i++){
				if(empirical.get(i)>0 || sampled.get(i)>0){hist.print(b.clear().append(i).tab().append(empirical.get(i)).tab().append(sampled.get(i)).nl());}
			}
			summary.println("arm\teligible_genomes\tsampled_genomes\tshreds\tcore_bases\tcontext_bases\tseed\tflank\tpolicy\tlabel");
			summary.print(b.clear().append(arm).tab().append(eligibleGenomes).tab().append(genomeCount).tab().append(written).tab().append(coreBases).tab().append(contextBases).tab().append(seed).tab().append(flank)
				.append("\tBALANCED_GENOMES_EMPIRICAL_LENGTH_UNIFORM_ELIGIBLE_START_WITH_REPLACEMENT\tSOURCE_MASKED_CM_SCREEN_PENDING\n"));
		}finally{final boolean failed=core.close()|context.close();RrnaResourceIO.finish(rows, genomes, hist, summary);require(!failed, "Background FASTA write failed");}
		System.err.println("COBALAMIN_SHRED_SAMPLE_PASS arm="+arm+" shreds="+written+" bases="+coreBases+"; all candidates still require both-strand clan CM screening");
	}
	private static ArrayList<Start> startSpace(ArrayList<Contig> contigs, int length){
		require(length>0, "Complement start-space requires positive fixed core length");final ArrayList<Start> result=new ArrayList<Start>();long total=0;
		for(Contig c:contigs){
			final long last=(long)c.bases.length-length+1;if(last<1){continue;}long next=1;
			for(int i=0; i<c.mask.starts.size; i++){
				final long from=Math.max(1L, (long)c.mask.starts.get(i)-length+1), to=Math.min(last, c.mask.ends.get(i));
				if(from>next){final long end=Math.min(last, from-1);if(end>=next){total=Math.addExact(total, end-next+1);result.add(new Start(c, (int)next, total));}}
				next=Math.max(next, to+1);if(next>last){break;}
			}
			if(next<=last){total=Math.addExact(total, last-next+1);result.add(new Start(c, (int)next, total));}
		}return result;
	}
	private static String stem(String assembly){require(!assembly.isEmpty(), "Assembly provenance cannot be empty");final int dot=assembly.indexOf('.');return dot<0?assembly:assembly.substring(0, dot);}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	private static final class Genome {
		Genome(String n, String a, String t, String f){name=n;assembly=a;taxid=t;fasta=f;}
		final String name, assembly, taxid, fasta;boolean panel;int sourceLoci;
		final ObjectMap<String,Mask> masks=new ObjectMap<String,Mask>(String.class, Mask.class);final ArrayList<String> maskNames=new ArrayList<String>();
	}
	private static final class Mask {Mask(int n){length=n;}final int length;final ArrayList<int[]> spans=new ArrayList<int[]>();}
	private static final class Contig {Contig(String n, byte[] b, Euk18sShredSet.Contig m){name=n;bases=b;mask=m;}final String name;final byte[] bases;final Euk18sShredSet.Contig mask;}
	private static final class Start {Start(Contig c, int s, long e){contig=c;start=s;end=e;}final Contig contig;final int start;final long end;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	private final ArrayList<Genome> selected=new ArrayList<Genome>();private final ObjectMap<String,Genome> byGenome=new ObjectMap<String,Genome>(String.class, Genome.class);
	private final IntList lengths=new IntList();private int[] schedule;private String arm;
	private int count=-1, flank=-1, genomeCount=-1, expectedGenomes=-1, expectedPanel=-1, expectedLengths=-1, expectedLoci=-1, eligibleGenomes;
	private long seed=Long.MIN_VALUE;
	private static final String[] FAMILIES={"RF00174", "RF01482", "RF01689"};
	public static final String HEADER="shred\tarm\tgenome\tassembly\ttaxid\tseqid\tcore_start1\tcore_stop1\tcontext_start1\tcontext_stop1\tcore_in_context_start1\tcore_in_context_stop1\tlength";
}
