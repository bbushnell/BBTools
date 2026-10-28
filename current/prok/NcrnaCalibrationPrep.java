package prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import shared.Tools;
import stream.Read;
import structures.ByteBuilder;

/** Freezes coordinate views and actual caller windows from an attributed corpus.
 * Source annotations, full-model CM cores and unresolved regions stay separate.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaCalibrationPrep {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args, NcrnaCalibrationPrep.class, false);
		String roster=null, audit=null, stages=null, out=null, family=null;
		int genomes=-1, loci=-1, windows=-1, contigs=-1, clen=-1;double ga=Double.NaN;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0, "Use flag=value");
			final String key=arg.substring(0, eq).toLowerCase(Locale.ROOT), value=arg.substring(eq+1);
			if(key.equals("roster")){roster=value;}else if(key.equals("audit")){audit=value;}
			else if(key.equals("stages")){stages=value;}else if(key.equals("outprefix")){out=value;}
			else if(key.equals("family")){family=value;}else if(key.equals("genomes")){genomes=Integer.parseInt(value);}
			else if(key.equals("loci")){loci=Integer.parseInt(value);}else if(key.equals("windows")){windows=Integer.parseInt(value);}
			else if(key.equals("contigs")){contigs=Integer.parseInt(value);}else if(key.equals("clen")){clen=Integer.parseInt(value);}
			else if(key.equals("gabits")){ga=Double.parseDouble(value);}else{throw new IllegalArgumentException("Unknown option: "+key);}
		}
		require(roster!=null && audit!=null && stages!=null && out!=null && family!=null && genomes>0 && loci>0
			&& windows>=0 && contigs>0 && clen>0 && Double.isFinite(ga), "Bind the complete cohort, audit, caller stages, model criteria and denominators");
		Shared.setThreads(1);Shared.TRIM_READ_DESCRIPTION=false;
		final NcrnaCalibrationPrep tool=new NcrnaCalibrationPrep(family);
		tool.loadRoster(roster, genomes);tool.loadAudit(audit, loci, clen, ga);tool.loadStages(stages, windows);
		tool.write(out, contigs);Shared.closeStream(pp.outstream);
	}

	private NcrnaCalibrationPrep(final String family_){family=family_;}

	/*--------------------------------------------------------------*/
	/*----------------            Inputs            ----------------*/
	/*--------------------------------------------------------------*/

	private void loadRoster(final String path, final int expected){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{
			RrnaResourceIO.header(in, ROSTER_HEADER, path);
			for(byte[] row; (row=in.nextLine())!=null; ){
				p.set(row);require(p.terms()==11, "Genome roster schema changed");
				final Genome g=new Genome(p.parseString(1), p.parseString(7));
				require(!genomeIds.contains(g.id) && !g.id.isEmpty() && !g.fasta.isEmpty(), "Every cohort genome needs one explicit FASTA: "+g.id);
				genomeIds.put(g.id, g);genomes.add(g);
			}
		}finally{require(!in.close(), "Genome roster read failed");}
		require(genomes.size()==expected, "Genome denominator changed");
	}
	private void loadAudit(final String path, final int expected, final int clen, final double ga){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		final ObjectSet<String> ids=new ObjectSet<String>(String.class);
		try{
			RrnaResourceIO.header(in, NcrnaCoreRecovery.AUDIT_HEADER, path);
			for(byte[] row; (row=in.nextLine())!=null; ){
				p.set(row);require(p.terms()==34, "Expected complete source plus CM-attribution fields");
				final String id=p.parseString(0), genome=p.parseString(1), seqid=p.parseString(6), strand=p.parseString(7), interval=p.parseString(8);
				require(ids.add(id) && genomeIds.contains(genome) && (strand.equals("+") || strand.equals("-")), "Duplicate locus, foreign genome or invalid strand: "+id);
				final int dash=interval.indexOf('-');require(dash>0 && interval.indexOf(',')<0, "Calibration cohort requires contiguous source loci: "+id);
				final int start=Integer.parseInt(interval.substring(0, dash)), stop=Integer.parseInt(interval.substring(dash+1)), length=p.parseInt(15);
				require(start>0 && stop>=start && stop<=length && stop-start+1==p.parseInt(9), "Source coordinates must match the declared contig and length: "+id);
				final Contig c=contig(seqid, length);
				require(c.genome==null || c.genome.equals(genome), "A source contig ID cannot ambiguously identify two genomes: "+seqid);c.genome=genome;
				final Locus l=new Locus(id, c, strand.charAt(0), start, stop, p.parseString(21));
				if(l.state.equals("CM_FULL_SPAN_SUPPORTED")){
					final int left=p.parseInt(18), right=p.parseInt(19), from=p.parseInt(26), to=p.parseInt(27);
					require(left>=0 && right>=0 && from>0 && to>=from && to<=stop-start+1+left+right
						&& p.parseInt(23)==1 && p.parseInt(24)==1 && p.parseInt(25)==clen && p.termEquals("no", 28) && p.parseDouble(29)>=ga,
						"Full-model view requires one overlapping, untruncated, complete, GA-passing hit: "+id);
					require(p.termEquals("EXACT", 33) || p.termEquals("HIT_WITHIN_ANNOTATION", 33) || p.termEquals("ANNOTATION_WITHIN_HIT", 33),
						"Full-supported coordinate views require containment, not an unrelated neighboring hit: "+id);
					l.coreStart=l.strand=='+'?start-left+from-1:stop+left-to+1;
					l.coreStop=l.strand=='+'?start-left+to-1:stop+left-from+1;
					require(l.coreStart>0 && l.coreStop<=length && l.coreStop-l.coreStart==to-from, "CM genomic projection must preserve length and fit its contig: "+id);
				}
				loci.add(l);
			}
		}finally{require(!in.close(), "CM audit read failed");}
		require(loci.size()==expected, "Locus denominator changed");
	}
	private void loadStages(final String path, final int expected){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{
			RrnaResourceIO.header(in, NcrnaStageDiagLogSink.HEADER, path);
			for(byte[] row; (row=in.nextLine())!=null; ){
				p.set(row);require(p.terms()==21, "Caller stage schema changed");
				if(!p.termEquals(family, 1) || !p.termEquals("WINDOW", 0)){continue;}
				final String seqid=RrnaResourceIO.first(p.parseString(2));final int length=p.parseInt(4);
				final Contig c=contig(seqid, length);
				final Window w=new Window(windows.size(), c, p.parseInt(3), p.parseInt(5), p.parseInt(6), p.parseInt(7),
					p.parseString(8), p.parseInt(10), p.parseInt(11), p.parseFloat(13));
				require((w.strand==0 || w.strand==1) && (w.pass==1 || w.pass==2) && w.start>=0 && w.stop>=w.start && w.stop<length,
					"Caller window must remain in its actual oriented contig");
				c.windows.add(w);windows.add(w);
			}
		}finally{require(!in.close(), "Caller stage read failed");}
		require(windows.size()==expected, "Window denominator changed: "+windows.size()+" expected="+expected);
	}
	private Contig contig(final String name, final int length){
		require(length>0, "Contig length must be positive");
		Contig c=contigs.get(name);
		if(c==null){c=new Contig(name, length);contigs.put(name, c);}
		require(c.length==length, "Conflicting lengths for contig: "+name);return c;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Outputs           ----------------*/
	/*--------------------------------------------------------------*/

	private void write(final String prefix, final int expectedContigs){
		RrnaResourceIO.fresh(prefix, ".windows.fa.gz", ".windows.tsv", ".source.tsv", ".full.tsv", ".unresolved.tsv", ".loci.tsv", ".summary.tsv");
		final ByteStreamWriter fasta=RrnaResourceIO.writer(prefix+".windows.fa.gz");final ByteBuilder b=new ByteBuilder();
		final ObjectSet<String> allIds=new ObjectSet<String>(String.class);final int[] emitted={0};
		try{
			for(Genome g:genomes){
				RrnaResourceIO.stream(g.fasta, r->{
					final String name=RrnaResourceIO.first(r.id);
					require(allIds.add(name), "Combined GFF grading requires globally unique cohort contig IDs: "+name);
					final Contig c=contigs.get(name);if(c==null){return;}
					require(!c.seen && r.bases!=null && r.length()==c.length && (c.genome==null || c.genome.equals(g.id)),
						"Stage/source contig must bind exactly once to its original genome and length: "+name);
					c.seen=true;c.genome=g.id;Tools.toUpperCase(r.bases);
					for(Window w:c.windows){
						final int first=w.strand==0?w.start:c.length-w.stop-1, last=w.strand==0?w.stop:c.length-w.start-1;
						final byte[] query=Arrays.copyOfRange(r.bases, first, last+1);
						if(w.strand==1){AminoAcid.reverseComplementBasesInPlace(query);}
						fasta.print(b.clear().append(">w").append(w.id).nl().append(query).nl());emitted[0]++;
					}
				});
			}
		}finally{RrnaResourceIO.finish(fasta);}
		require(allIds.size()==expectedContigs && emitted[0]==windows.size(), "Genome/window extraction must conserve every frozen contig and window");
		for(String key:contigs.keys()){if(key!=null){require(contigs.get(key).seen, "Source/stage contig is absent from the frozen genomes: "+key);}}
		final ByteStreamWriter source=RrnaResourceIO.writer(prefix+".source.tsv"), full=RrnaResourceIO.writer(prefix+".full.tsv"), unresolved=RrnaResourceIO.writer(prefix+".unresolved.tsv");
		final ByteStreamWriter ledger=RrnaResourceIO.writer(prefix+".loci.tsv"), rows=RrnaResourceIO.writer(prefix+".windows.tsv");int supported=0;
		try{
			source.println(TRUTH_HEADER);full.println(TRUTH_HEADER);unresolved.println(TRUTH_HEADER);
			ledger.println("locus_id\tgenome\tseqid\tstrand\tsource_start1\tsource_stop1\tcm_start1\tcm_stop1\tcm_state");
			for(Locus l:loci){
				truth(source, b, l, l.start, l.stop);
				if(l.coreStart>0){truth(full, b, l, l.coreStart, l.coreStop);supported++;}else{truth(unresolved, b, l, l.start, l.stop);}
				ledger.print(b.clear().append(l.id).tab().append(l.contig.genome).tab().append(l.contig.name).tab().append(l.strand).tab()
					.append(l.start).tab().append(l.stop).tab().append(l.coreStart).tab().append(l.coreStop).tab().append(l.state).nl());
			}
			rows.println(WINDOW_HEADER);
			for(Window w:windows){
				rows.print(b.clear().append('w').append(w.id).tab().append(w.contig.genome).tab().append(w.contig.name).tab().append(w.strand)
					.tab().append(w.contig.length).tab().append(w.pass).tab().append(w.start).tab().append(w.stop).tab().append(w.outcome)
					.tab().append(w.shortlisted).tab().append(w.model).tab().appendSlow(w.identity).tab().append(w.start==0 || w.stop==w.contig.length-1).nl());
			}
		}finally{RrnaResourceIO.finish(source, full, unresolved, ledger, rows);}
		final ByteStreamWriter summary=RrnaResourceIO.writer(prefix+".summary.tsv");
		try{
			summary.println("genomes\tcontigs\tsource_loci\tfull_model_loci\tunresolved_loci\twindows");
			summary.print(b.clear().append(genomes.size()).tab().append(allIds.size()).tab().append(loci.size()).tab().append(supported).tab()
				.append(loci.size()-supported).tab().append(windows.size()).nl());
		}finally{RrnaResourceIO.finish(summary);}
		System.err.println("NCRNA_CALIBRATION_PREP_PASS genomes="+genomes.size()+" source="+loci.size()+" full="+supported+" windows="+windows.size());
	}
	private void truth(final ByteStreamWriter out, final ByteBuilder b, final Locus l, final int start, final int stop){
		assert(start>0 && stop>=start):"Truth views use forward one-based inclusive coordinates";
		out.print(b.clear().append(l.contig.name).tab().append(family).tab().append(start).tab().append(stop).tab().append(l.strand).nl());
	}
	private static void require(final boolean ok, final String why){if(!ok){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Genome {
		Genome(String id_, String fasta_){id=id_;fasta=fasta_;}
		final String id, fasta;
	}
	private static final class Contig {
		Contig(String name_, int length_){name=name_;length=length_;}
		final String name;final int length;String genome;boolean seen;
		final ArrayList<Window> windows=new ArrayList<Window>();
	}
	private static final class Locus {
		Locus(String id_, Contig c, char strand_, int start_, int stop_, String state_){id=id_;contig=c;strand=strand_;start=start_;stop=stop_;state=state_;}
		final String id, state;final Contig contig;final char strand;final int start, stop;int coreStart, coreStop;
	}
	private static final class Window {
		Window(int id_, Contig c, int strand_, int pass_, int start_, int stop_, String outcome_, int shortlisted_, int model_, float identity_){
			id=id_;contig=c;strand=strand_;pass=pass_;start=start_;stop=stop_;outcome=outcome_;shortlisted=shortlisted_;model=model_;identity=identity_;
		}
		final Contig contig;final String outcome;final int id, strand, pass, start, stop, shortlisted, model;final float identity;
	}

	private final String family;
	private final ArrayList<Genome> genomes=new ArrayList<Genome>();
	private final ObjectMap<String,Genome> genomeIds=new ObjectMap<String,Genome>(String.class, Genome.class);
	private final ObjectMap<String,Contig> contigs=new ObjectMap<String,Contig>(String.class, Contig.class);
	private final ArrayList<Locus> loci=new ArrayList<Locus>();
	private final ArrayList<Window> windows=new ArrayList<Window>();
	public static final String ROSTER_HEADER="index\tgenome\tassembly\ttaxid\tdomain\tphylum\tgff\tfasta\tloci\tgenome_bases\tcontigs";
	public static final String TRUTH_HEADER="seqid\tfamily\tstart\tstop\tstrand";
	public static final String WINDOW_HEADER="window\tgenome\tseqid\tstrand\tcontig_length\tpass\tstart0\tstop0\toutcome\tshortlist_size\tmodel\tidentity\tphysical_edge";
}
