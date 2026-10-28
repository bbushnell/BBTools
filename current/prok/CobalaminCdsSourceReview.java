package prok;

import java.util.ArrayList;
import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/** Binds removed CDS rows to original annotations by strand and translation stop.
 * The translation stop is genomic stop on +, genomic start on -. Annotated start
 * changes remain distinct from whole-gene loss; source presence is not function.
 * The producer verifies the supplied original GFF hashes before execution.
 * @author Brian Bushnell, Raiden
 */
public final class CobalaminCdsSourceReview {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args, CobalaminCdsSourceReview.class, false);
		String changes=null, contigs=null, roster=null, hashes=null, out=null;int removed=-1, added=-1;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0, "Use flag=value");
			final String key=arg.substring(0, eq).toLowerCase(Locale.ROOT), value=arg.substring(eq+1);
			if(key.equals("changes")){changes=value;}else if(key.equals("contigs")){contigs=value;}
			else if(key.equals("roster")){roster=value;}else if(key.equals("gffhashes")){hashes=value;}
			else if(key.equals("outprefix")){out=value;}else if(key.equals("removed")){removed=Integer.parseInt(value);}
			else if(key.equals("added")){added=Integer.parseInt(value);}else{throw new IllegalArgumentException(key);}
		}
		require(changes!=null && contigs!=null && roster!=null && hashes!=null && out!=null && removed>=0 && added>=0,
			"Bind the native gene delta, frozen panel/roster and verified source hashes with explicit counts");
		Shared.setThreads(1);
		try{
			final CobalaminCdsSourceReview x=new CobalaminCdsSourceReview();x.contigs(contigs);x.changes(changes, removed, added);
			x.roster(roster, hashes);x.annotations();x.write(out);
		}finally{Shared.closeStream(pp.outstream);}
	}
	private void contigs(String path){
		read(path, CobalaminPanelPrep.CONTIG_HEADER, p->{
			final String seq=p.parseString(2), genome=p.parseString(0);final int length=p.parseInt(3);
			require(!seq.isEmpty() && !genome.isEmpty() && length>0 && !byContig.contains(seq), "Contig ownership must remain unique");
			byContig.put(seq, new Contig(genome, length));
		});
	}
	private void changes(String path, int removedCount, int addedCount){
		final long[] line={1};
		final ObjectSet<String> anchors=new ObjectSet<String>(String.class);
		read(path, GffCallSetDelta.CHANGE_HEADER, p->{
			line[0]++;require(p.termEquals("REMOVED", 0) || p.termEquals("ADDED", 0), "Unknown native delta direction");
			if(!p.termEquals("CDS", 3)){return;}
			final Feature f=new Feature(p, line[0]);final Contig c=byContig.get(f.seqid);
			require(c!=null && f.stop<=c.length, "Changed CDS must bind the original frozen contig");f.genome=c.genome;
			if(p.termEquals("ADDED", 0)){added.add(f);}else{
				require(anchors.add(f.seqid+'\t'+f.plus+'\t'+f.translationStop()),
					"Multiple removed CDS rows at one translation stop require explicit event disambiguation");
				final Gene gene=new Gene(f);removed.add(gene);c.genes.add(gene);
				if(!byGenome.contains(c.genome)){final Genome g=new Genome(c.genome);byGenome.put(c.genome, g);genomes.add(g);}
			}
		});
		require(removed.size()==removedCount && added.size()==addedCount, "Every removed/added CDS row must survive the source review");
		for(Gene g:removed){for(Feature f:added){if(g.old.sameStop(f)){
			require(g.replacement==null && !f.paired, "Ambiguous stop-anchored CDS replacements require explicit review, not an arbitrary pairing");
			g.replacement=f;f.paired=true;
		}}}
	}
	private void roster(String path, String hashes){
		final ObjectMap<String,String> digest=EukPanelGffValidation.hashes(hashes);final ObjectSet<String> seen=new ObjectSet<String>(String.class);
		read(path, NcrnaCalibrationPrep.ROSTER_HEADER, p->{
			final String name=p.parseString(1);require(seen.add(name), "Original genome roster contains a repeated identity");
			final Genome g=byGenome.get(name);if(g==null){return;}
			// NcrnaCalibrationPrep.ROSTER_HEADER orders gff before fasta.
			g.gff=p.parseString(6);g.hash=digest.get(g.gff);
			require(g.hash!=null && g.hash.matches("[0-9a-f]{20}"), "Changed CDS genome needs an explicitly verified original GFF sha80");
		});
		for(Genome g:genomes){require(g.gff!=null, "Changed CDS lacks its original source genome: "+g.name);}
	}
	private void annotations(){
		for(Genome g:genomes){
			final ByteFile in=RrnaResourceIO.open(g.gff);final LineParser1 p=new LineParser1('\t');String last=null;Contig c=null;
			try{for(byte[] row; (row=in.nextLine())!=null; ){
				g.lines++;if(row.length==0){continue;}if(row[0]=='#'){if(Tools.startsWith(row, "##FASTA")){break;}continue;}
				p.set(row);require(p.terms()==9, "Original source GFF9 schema changed: "+g.gff);
				if(last==null || !p.termEquals(last, 0)){last=p.parseString(0);c=byContig.get(last);}
				require(c!=null && c.genome.equals(g.name), "Original GFF row belongs to a foreign or missing panel contig");
				if(!p.termEquals("CDS", 2)){continue;}g.cds++;
				if(c.genes.isEmpty()){continue;}
				final int start=p.parseInt(3), stop=p.parseInt(4);final boolean plus=p.termEquals("+", 6), minus=p.termEquals("-", 6);
				require(start>0 && start<=c.length && stop>=start && (long)stop-start+1<=c.length, "Original CDS spans at most one physical contig turn");
				if(!plus && !minus){continue;}
				for(Gene gene:c.genes){
					if(gene.old.plus==plus && gene.old.translationStop()==(plus?stop:start)){
						require(stop<=c.length, "A matched origin-crossing CDS needs explicit circular start-error handling before classification");
						gene.annotations.add(new Annotation(start, stop, g.lines, row.clone()));
					}
				}
			}}finally{require(!in.close(), "Source GFF read failed: "+g.gff);}
		}
	}
	private void write(String out){
		RrnaResourceIO.fresh(out, ".summary.tsv", ".annotations.tsv", ".added.tsv", ".sources.tsv");
		final ByteStreamWriter summary=RrnaResourceIO.writer(out+".summary.tsv"), detail=RrnaResourceIO.writer(out+".annotations.tsv");
		final ByteStreamWriter additions=RrnaResourceIO.writer(out+".added.tsv"), sources=RrnaResourceIO.writer(out+".sources.tsv");final ByteBuilder b=new ByteBuilder();
		try{
			summary.println("delta_line\tgenome\tseqid\tstrand\told_start1\told_stop1\tnew_start1\tnew_stop1\tchange\tsource_CDS_stop_matches\tclassification\tstart_improved_matches\tstart_worsened_matches\tstart_equal_error_matches");
			detail.println("delta_line\tgenome\tgff\tgff_sha80\tannotation_line\ttranslation_stop1\told_start_offset_bp\tnew_start_offset_bp\tcomparison\tseqid\tsource\ttype\tstart1\tstop1\tscore\tstrand\tphase\tattributes");
			additions.println("delta_line\tgenome\tseqid\tstart1\tstop1\tstrand\tpaired_same_translation_stop");
			sources.println("genome\tgff\tgff_sha80\tlines_read\tCDS_rows_read");
			for(Gene gene:removed){
				final Feature old=gene.old, replacement=gene.replacement;final Genome genome=byGenome.get(old.genome);
				int improved=0, worsened=0, equal=0;
				for(Annotation a:gene.annotations){
					final int oldOffset=old.plus?old.start-a.start:a.stop-old.stop;
					final int newOffset=replacement==null?0:old.plus?replacement.start-a.start:a.stop-replacement.stop;
					final String comparison;
					if(replacement==null){comparison="NO_REPLACEMENT";}
					else if(Math.abs(newOffset)<Math.abs(oldOffset)){improved++;comparison="IMPROVED";}
					else if(Math.abs(newOffset)>Math.abs(oldOffset)){worsened++;comparison="WORSENED";}
					else{equal++;comparison="SAME_ERROR";}
					detail.print(b.clear().append(old.line).tab().append(old.genome).tab().append(genome.gff).tab().append(genome.hash).tab().append(a.line)
						.tab().append(old.translationStop()).tab().append(oldOffset).tab().append(replacement==null?"NA":Integer.toString(newOffset))
						.tab().append(comparison).tab().append(a.raw).nl());
				}
				final String classification=gene.annotations.isEmpty()?"PREDICTED-ONLY":replacement==null?"REAL-GENE-LOSS":"ANNOTATED-START-CHANGE";
				summary.print(b.clear().append(old.line).tab().append(old.genome).tab().append(old.seqid).tab().append(old.plus?'+':'-')
					.tab().append(old.start).tab().append(old.stop).tab().append(replacement==null?"NA":Integer.toString(replacement.start))
					.tab().append(replacement==null?"NA":Integer.toString(replacement.stop)).tab().append(replacement==null?"REMOVED":"SAME_STOP_START_CHANGED")
					.tab().append(gene.annotations.size()).tab().append(classification).tab().append(improved).tab().append(worsened).tab().append(equal).nl());
			}
			for(Feature f:added){additions.print(b.clear().append(f.line).tab().append(f.genome).tab().append(f.seqid).tab().append(f.start).tab().append(f.stop).tab().append(f.plus?'+':'-').tab().append(f.paired).nl());}
			for(Genome g:genomes){sources.print(b.clear().append(g.name).tab().append(g.gff).tab().append(g.hash).tab().append(g.lines).tab().append(g.cds).nl());}
		}finally{RrnaResourceIO.finish(summary, detail, additions, sources);}
		System.err.println("COBALAMIN_CDS_SOURCE_REVIEW_PASS removed="+removed.size()+" added="+added.size()+"; same translation-stop/strand annotation rule, not an independent functional assay");
	}
	private static void read(String path, String header, java.util.function.Consumer<LineParser1> action){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');p.set(header.getBytes(java.nio.charset.StandardCharsets.US_ASCII));final int columns=p.terms();
		try{RrnaResourceIO.header(in, header, path);for(byte[] row; (row=in.nextLine())!=null; ){p.set(row);require(p.terms()==columns, "Native table schema changed: "+path);action.accept(p);}}
		finally{require(!in.close(), "Native table read failed: "+path);}
	}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final class Feature {
		Feature(LineParser1 p, long n){
			line=n;seqid=p.parseString(1);start=p.parseInt(4);stop=p.parseInt(5);plus=p.termEquals("+", 7);
			require(start>0 && stop>=start && (plus || p.termEquals("-", 7)), "Changed CDS needs positive coordinates and a known strand");
		}
		int translationStop(){return plus?stop:start;}
		boolean sameStop(Feature other){return seqid.equals(other.seqid) && plus==other.plus && translationStop()==other.translationStop();}
		final long line;final String seqid;final int start, stop;final boolean plus;String genome;boolean paired;
	}
	private static final class Gene {Gene(Feature f){old=f;}final Feature old;Feature replacement;final ArrayList<Annotation> annotations=new ArrayList<Annotation>();}
	private static final class Annotation {Annotation(int a, int z, long n, byte[] r){start=a;stop=z;line=n;raw=r;}final int start, stop;final long line;final byte[] raw;}
	private static final class Contig {Contig(String g, int n){genome=g;length=n;}final String genome;final int length;final ArrayList<Gene> genes=new ArrayList<Gene>();}
	private static final class Genome {Genome(String n){name=n;}final String name;String gff, hash;long lines, cds;}
	private final ObjectMap<String,Contig> byContig=new ObjectMap<String,Contig>(String.class, Contig.class);
	private final ObjectMap<String,Genome> byGenome=new ObjectMap<String,Genome>(String.class, Genome.class);
	private final ArrayList<Genome> genomes=new ArrayList<Genome>();
	private final ArrayList<Gene> removed=new ArrayList<Gene>();private final ArrayList<Feature> added=new ArrayList<Feature>();
}
