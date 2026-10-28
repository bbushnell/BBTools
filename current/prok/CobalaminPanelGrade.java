package prok;

import java.util.ArrayList;
import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.PreParser;
import prok.NcrnaGffFamilyGrader.Locus;
import prok.NcrnaGffFamilyGrader.MatchResult;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;
import structures.IntList;

/** One-to-one B12 subtype/union grades with explicit unresolved subtype abstentions.
 * Historical coordinates and reviewed full-core coordinates are separate views.
 * @author Brian Bushnell, Raiden
 */
public final class CobalaminPanelGrade {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args, CobalaminPanelGrade.class, false);
		String panel=null, contigs=null, calls=null, out=null;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0, "Use flag=value");
			final String k=arg.substring(0, eq).toLowerCase(Locale.ROOT), v=arg.substring(eq+1);
			if(k.equals("panel")){panel=v;}else if(k.equals("contigs")){contigs=v;}
			else if(k.equals("in")){calls=v;}else if(k.equals("outprefix")){out=v;}else{throw new IllegalArgumentException("Unknown option: "+k);}
		}
		require(panel!=null && contigs!=null && calls!=null && out!=null, "Bind reviewed panel, exact contig ownership and actual caller GFF");
		Shared.setThreads(1);
		try{final CobalaminPanelGrade x=new CobalaminPanelGrade();x.contigs(contigs);x.truth(panel);x.calls(calls);x.write(out);}
		finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Inputs            ----------------*/
	/*--------------------------------------------------------------*/

	private void contigs(String path){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{
			RrnaResourceIO.header(in, CobalaminPanelPrep.CONTIG_HEADER, path);
			for(byte[] row; (row=in.nextLine())!=null; ){
				p.set(row);require(p.terms()==4, "Panel contig ownership schema changed");
				final String id=p.parseString(0), assembly=p.parseString(1), name=p.parseString(2);final int length=p.parseInt(3);
				Genome g=genomeMap.get(id);if(g==null){g=new Genome(id, assembly, genomes.size());genomeMap.put(id, g);genomes.add(g);}
				require(g.assembly.equals(assembly) && length>0 && !contigMap.contains(name), "Contig IDs must identify one exact genome and length");
				contigMap.put(name, new Contig(g, name, length));
			}
		}finally{require(!in.close(), "Contig ownership read failed");}
		require(!genomes.isEmpty(), "A grade requires its complete genome denominator");
	}
	private void truth(String path){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		final ObjectSet<String> seen=new ObjectSet<String>(String.class), physical=new ObjectSet<String>(String.class), cores=new ObjectSet<String>(String.class);
		try{
			RrnaResourceIO.header(in, CobalaminPanelPrep.HEADER, path);
			for(byte[] row; (row=in.nextLine())!=null; ){
				p.set(row);require(p.terms()==10, "Reviewed panel schema changed");
				final String id=p.parseString(0), family=p.parseString(3), member=p.parseString(4);final int source=family(family);
				final Contig c=contigMap.get(p.parseString(2));
				require(c!=null && p.termEquals(c.genome.id, 1) && source>=0 && seen.add(family+'/'+id), "Truth must retain unique IDs and exact genome ownership");
				final Locus original=locus(c, family, p.parseInt(5), p.parseInt(6), p.parseString(9));
				require(physical.add(key(original)), "Exact duplicate source coordinates require independent adjudication");
				final int start=p.parseInt(7), stop=p.parseInt(8);final Locus core;
				if(member.equals(CobalaminPanelPrep.UNRESOLVED)){require(start==0 && stop==0, "Unresolved truth cannot carry a claimed full core");core=null;}
				else{
					require(member.equals(family) || member.equals(CobalaminPanelPrep.CLAN_ONLY), "Only reviewed own-subtype or clan-only labels are accepted");
					core=locus(c, member, start, stop, p.parseString(9));
					require(cores.add(key(core)), "Distinct source intervals share an exact reviewed core; adjudication is required before FULL_CM union grading: "+id);
				}
				truth.add(new Truth(id, c, source, original, core, member.equals(CobalaminPanelPrep.CLAN_ONLY)));
			}
		}finally{require(!in.close(), "Reviewed panel read failed");}
	}
	private void calls(String path){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		final ObjectSet<String> seen=new ObjectSet<String>(String.class);
		try{for(byte[] row; (row=in.nextLine())!=null; ){
			if(row.length==0 || row[0]=='#'){continue;}p.set(row);require(p.terms()==9, "Malformed caller GFF");
			p.setBounds(8);if(Tools.indexOf(row, "cobalamin", p.a())<0){continue;}
			final String attrs=p.parseString(8), model=attribute(attrs, "model"), name=attribute(attrs, "family");
			final boolean selected=model!=null && model.startsWith("cobalamin") || name!=null && name.startsWith("cobalamin");
			if(!selected){continue;}
			require(p.termEquals("RNA", 2) && model!=null && name!=null, "Every B12 call requires RNA type, model and family attribution");
			int index=-1;for(int i=0; i<FAMILIES.length; i++){
				if(name.equals("cobalamin_"+FAMILIES[i].toLowerCase(Locale.ROOT))){index=i;break;}
			}
			require(index>=0 && model.startsWith("cobalamin_"+FAMILIES[index]+'_'), "Unknown or conflicting B12 subtype attribution: "+attrs);
			final Contig c=contigMap.get(p.parseString(0));require(c!=null, "Call contig is outside the frozen panel");
			final Locus l=locus(c, FAMILIES[index], p.parseInt(3), p.parseInt(4), p.parseString(6));
			calls.add(new Call(l, c, index, model, !seen.add(key(l)+'\t'+index)));
		}}finally{require(!in.close(), "Caller GFF read failed");}
	}
	private static String attribute(String attrs, String name){
		final String needle=name+':';String value=null;
		for(int pos=attrs.indexOf(needle); pos>=0; pos=attrs.indexOf(needle, pos+needle.length())){
			if(pos>0 && attrs.charAt(pos-1)!=',' && attrs.charAt(pos-1)!=';'){continue;}
			require(value==null, "Duplicate caller attribute: "+name);int end=pos+needle.length();
			while(end<attrs.length() && attrs.charAt(end)!=',' && attrs.charAt(end)!=';'){end++;}
			value=attrs.substring(pos+needle.length(), end);require(!value.isEmpty(), "Empty caller attribute: "+name);
		}return value;
	}
	private static Locus locus(Contig c, String family, int start, int stop, String strand){
		require((strand.equals("+") || strand.equals("-")) && start>0 && stop>=start && stop<=c.length, "Locus must remain on its exact contig and strand");
		return new Locus(c.id, family, start, stop, strand.charAt(0));
	}
	private static String key(Locus l){return l.seqid+'\t'+l.strand+'\t'+l.start+'\t'+l.stop;}
	private static int family(String name){for(int i=0; i<FAMILIES.length; i++){if(FAMILIES[i].equals(name)){return i;}}return -1;}

	/*--------------------------------------------------------------*/
	/*----------------           Grading            ----------------*/
	/*--------------------------------------------------------------*/

	private void write(String prefix){
		RrnaResourceIO.fresh(prefix, ".summary.tsv", ".genomes.tsv", ".calls.tsv", ".misses.tsv");
		final ByteStreamWriter summary=RrnaResourceIO.writer(prefix+".summary.tsv"), perGenome=RrnaResourceIO.writer(prefix+".genomes.tsv");
		final ByteStreamWriter audit=RrnaResourceIO.writer(prefix+".calls.tsv"), misses=RrnaResourceIO.writer(prefix+".misses.tsv");
		try{
			summary.println("view\tscope\tgenomes\ttruth\tcalls\tTP\tFP\tFN\tunscored\texact_duplicates\tprecision\trecall\tapparent_FP_per_genome\texact5\texact3\tboth_exact");
			perGenome.println("view\tscope\tgenome\ttruth\tcalls\tTP\tFP\tFN\tunscored\texact_duplicates\texact5\texact3\tboth_exact");
			audit.println("view\tscope\tcall_index\tgenome\tseqid\temitted_family\tstart1\tstop1\tstrand\tmodel\toutcome\ttruth_id\td5\td3");
			misses.println("view\tscope\ttruth_id\tgenome\tseqid\tstart1\tstop1\tstrand");
			for(boolean full:new boolean[]{false, true}){for(int scope=0; scope<4; scope++){grade(full, scope, summary, perGenome, audit, misses);}}
		}finally{RrnaResourceIO.finish(summary, perGenome, audit, misses);}
		System.err.println("COBALAMIN_PANEL_GRADE_PASS genomes="+genomes.size()+" truth="+truth.size()+" calls="+calls.size()+"; apparent metrics, unresolved evidence retained");
	}
	private void grade(boolean full, int scope, ByteStreamWriter summary, ByteStreamWriter perGenome, ByteStreamWriter audit, ByteStreamWriter misses){
		assert(scope>=0 && scope<=3):"Three biological subtypes and their pooled union are the only report scopes";
		final String view=full?"FULL_CM":"SOURCE", label=scope==3?"union":FAMILIES[scope];
		final ArrayList<Locus> refs=new ArrayList<Locus>(), scored=new ArrayList<Locus>(), masks=new ArrayList<Locus>();
		final ArrayList<Truth> owners=new ArrayList<Truth>();final IntList indices=new IntList();
		final Stats[] byGenome=new Stats[genomes.size()];for(int i=0; i<byGenome.length; i++){byGenome[i]=new Stats();}
		for(Truth t:truth){
			if(!full){if(scope==3 || t.source==scope){refs.add(t.original);owners.add(t);byGenome[t.contig.genome.index].truth++;}}
			else if(t.core==null || scope!=3 && t.clanOnly){masks.add(t.core==null?t.original:t.core);}
			else if(scope==3 || t.source==scope){refs.add(t.core);owners.add(t);byGenome[t.contig.genome.index].truth++;}
		}
		final boolean[] unscored=new boolean[calls.size()];final ObjectSet<String> unionSeen=new ObjectSet<String>(String.class);
		for(int i=0; i<calls.size(); i++){
			final Call c=calls.get(i);if(scope!=3 && c.family!=scope){continue;}
			final boolean duplicate=scope==3?!unionSeen.add(key(c.locus)):c.duplicate;
			final Stats s=byGenome[c.contig.genome.index];s.calls++;if(duplicate){s.duplicates++;}
			// A repeated prediction cannot consume a second overlapping reference.
			// Retain it for the audit/FP count below, but exclude it from matching.
			if(duplicate){continue;}
			boolean edge=false, masked=false;
			if(full && !duplicate){
				for(Locus r:refs){if(NcrnaGffFamilyGrader.overlapsReciprocal(c.locus, r, .5f)){edge=true;break;}}
				if(!edge){for(Locus mask:masks){if(overlap(c.locus, mask)){masked=true;break;}}}
			}
			if(masked){unscored[i]=true;s.unscored++;}else{scored.add(c.locus);indices.add(i);}
		}
		final MatchResult match=NcrnaGffFamilyGrader.matchOneToOne(scored, refs, .5f);final int[] truthForCall=new int[calls.size()];java.util.Arrays.fill(truthForCall, -1);
		for(int j=0; j<indices.size; j++){truthForCall[indices.get(j)]=match.callMatchedTruthIdx[j];}
		final ByteBuilder b=new ByteBuilder();
		for(int i=0; i<calls.size(); i++){
			final Call c=calls.get(i);if(scope!=3 && c.family!=scope){continue;}
			final Stats s=byGenome[c.contig.genome.index];final int ti=truthForCall[i];final String outcome=unscored[i]?"UNSCORED":ti>=0?"TP":"APPARENT_FP";
			if(!unscored[i]){if(ti>=0){s.tp++;}else{s.fp++;}}
			b.clear().append(view).tab().append(label).tab().append(i).tab().append(c.contig.genome.id).tab().append(c.locus.seqid).tab().append(FAMILIES[c.family])
				.tab().append(c.locus.start).tab().append(c.locus.stop).tab().append(c.locus.strand).tab().append(c.model).tab().append(outcome).tab();
			if(ti<0){b.append("NA\tNA\tNA");}else{
				final Locus r=refs.get(ti);final int d5=c.locus.strand=='+'?c.locus.start-r.start:r.stop-c.locus.stop;
				final int d3=c.locus.strand=='+'?c.locus.stop-r.stop:r.start-c.locus.start;
				if(d5==0){s.exact5++;}if(d3==0){s.exact3++;}if(d5==0 && d3==0){s.both++;}
				b.append(owners.get(ti).id).tab().append(d5).tab().append(d3);
			}audit.print(b.nl());
		}
		for(int i=0; i<refs.size(); i++){if(!match.truthConsumed[i]){
			final Truth t=owners.get(i);final Locus l=refs.get(i);byGenome[t.contig.genome.index].fn++;
			misses.print(b.clear().append(view).tab().append(label).tab().append(t.id).tab().append(t.contig.genome.id).tab().append(l.seqid).tab().append(l.start).tab().append(l.stop).tab().append(l.strand).nl());
		}}
		final Stats total=new Stats();
		for(int i=0; i<byGenome.length; i++){
			final Stats s=byGenome[i];s.check();total.add(s);
			perGenome.print(s.append(b.clear().append(view).tab().append(label).tab().append(genomes.get(i).id)).nl());
		}
		total.check();b.clear().append(view).tab().append(label).tab().append(genomes.size());
		b.tab().append(total.truth).tab().append(total.calls).tab().append(total.tp).tab().append(total.fp).tab().append(total.fn).tab().append(total.unscored).tab().append(total.duplicates);
		ratio(b, total.tp, total.tp+total.fp);ratio(b, total.tp, total.truth);ratio(b, total.fp, genomes.size());
		summary.print(b.tab().append(total.exact5).tab().append(total.exact3).tab().append(total.both).nl());
	}
	private static boolean overlap(Locus a, Locus b){return a.seqid.equals(b.seqid) && a.strand==b.strand && a.start<=b.stop && b.start<=a.stop;}
	private static void ratio(ByteBuilder b, long n, long d){b.tab();if(d==0){b.append("NA");}else{b.append(n/(double)d, 8);}}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final class Genome {
		Genome(String id_, String assembly_, int index_){id=id_;assembly=assembly_;index=index_;}final String id, assembly;final int index;
	}
	private static final class Contig {
		Contig(Genome g, String id_, int length_){genome=g;id=id_;length=length_;}final Genome genome;final String id;final int length;
	}
	private static final class Truth {
		Truth(String id_, Contig c, int source_, Locus original_, Locus core_, boolean clanOnly_){id=id_;contig=c;source=source_;original=original_;core=core_;clanOnly=clanOnly_;}
		final String id;final Contig contig;final int source;final Locus original, core;final boolean clanOnly;
	}
	private static final class Call {
		Call(Locus l, Contig c, int family_, String model_, boolean duplicate_){locus=l;contig=c;family=family_;model=model_;duplicate=duplicate_;}
		final Locus locus;final Contig contig;final int family;final String model;final boolean duplicate;
	}
	private static final class Stats {
		long truth, calls, tp, fp, fn, unscored, duplicates, exact5, exact3, both;
		void check(){require(truth==tp+fn && calls==tp+fp+unscored, "Truth/query conservation must hold per genome and globally");}
		void add(Stats s){truth+=s.truth;calls+=s.calls;tp+=s.tp;fp+=s.fp;fn+=s.fn;unscored+=s.unscored;duplicates+=s.duplicates;exact5+=s.exact5;exact3+=s.exact3;both+=s.both;}
		ByteBuilder append(ByteBuilder b){return b.tab().append(truth).tab().append(calls).tab().append(tp).tab().append(fp).tab().append(fn).tab().append(unscored).tab().append(duplicates).tab().append(exact5).tab().append(exact3).tab().append(both);}
	}
	private final ObjectMap<String,Genome> genomeMap=new ObjectMap<String,Genome>(String.class, Genome.class);
	private final ObjectMap<String,Contig> contigMap=new ObjectMap<String,Contig>(String.class, Contig.class);
	private final ArrayList<Genome> genomes=new ArrayList<Genome>();
	private final ArrayList<Truth> truth=new ArrayList<Truth>();private final ArrayList<Call> calls=new ArrayList<Call>();
	private static final String[] FAMILIES=CobalaminPanelPrep.FAMILIES;
}
