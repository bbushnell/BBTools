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
import structures.ByteBuilder;

/** Projects reviewed B12 evidence onto one exact genome panel without changing biological labels.
 * @author Brian Bushnell, Raiden
 */
public final class CobalaminPanelPrep {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args, CobalaminPanelPrep.class, false);
		String roster=null, legacy=null, rf14=null, rf16=null, out=null;
		int genomes=-1, contigs=-1, oldRows=-1, oldFull=-1, rows14=-1, rows16=-1;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0, "Use flag=value");
			final String k=arg.substring(0, eq).toLowerCase(Locale.ROOT), v=arg.substring(eq+1);
			if(k.equals("roster")){roster=v;}else if(k.equals("legacy")){legacy=v;}
			else if(k.equals("rf14")){rf14=v;}else if(k.equals("rf16")){rf16=v;}else if(k.equals("outprefix")){out=v;}
			else if(k.equals("genomes")){genomes=Integer.parseInt(v);}else if(k.equals("contigs")){contigs=Integer.parseInt(v);}
			else if(k.equals("legacyrows")){oldRows=Integer.parseInt(v);}else if(k.equals("legacyfull")){oldFull=Integer.parseInt(v);}
			else if(k.equals("rf14rows")){rows14=Integer.parseInt(v);}else if(k.equals("rf16rows")){rows16=Integer.parseInt(v);}
			else{throw new IllegalArgumentException("Unknown option: "+k);}
		}
		require(roster!=null && legacy!=null && rf14!=null && rf16!=null && out!=null && genomes>0 && contigs>0
			&& oldRows>=0 && oldFull>=0 && rows14>=0 && rows16>=0, "Bind the exact panel and complete reviewed evidence populations");
		Shared.setThreads(1);Shared.TRIM_READ_DESCRIPTION=false;
		try{
			final CobalaminPanelPrep x=new CobalaminPanelPrep();x.roster(roster, genomes, contigs);
			x.legacy(legacy, oldRows, oldFull);x.overlay(rf14, 1, rows14);x.overlay(rf16, 2, rows16);x.write(out);
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Inputs            ----------------*/
	/*--------------------------------------------------------------*/

	private void roster(String path, int expectedGenomes, int expectedContigs){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{
			RrnaResourceIO.header(in, NcrnaCalibrationPrep.ROSTER_HEADER, path);
			for(byte[] row; (row=in.nextLine())!=null; ){
				p.set(row);require(p.terms()==11, "Frozen panel roster schema changed");
				final Genome g=new Genome(p.parseString(1), p.parseString(2), p.parseString(7));
				require(!genomes.contains(g.id), "Duplicate panel genome: "+g.id);genomes.put(g.id, g);orderedGenomes.add(g);
				final long expectedBases=p.parseLong(9);final int expectedN=p.parseInt(10);
				RrnaResourceIO.stream(g.fasta, r->{
					final String name=RrnaResourceIO.first(r.id);
					require(r.bases!=null && r.length()>0 && !contigs.contains(name), "Combined grading requires unique nonempty contigs: "+name);
					final Contig c=new Contig(g, name, r.length());contigs.put(name, c);g.contigs.add(c);g.bases+=r.length();
				});
				// The caller-input roster uses zero/zero when original extraction skipped a target-free genome.
				require(expectedBases==0 && expectedN==0 || g.bases==expectedBases && g.contigs.size()==expectedN,
					"Measured source genome bases/contigs changed: "+g.id);
			}
		}finally{require(!in.close(), "Panel roster read failed");}
		require(genomes.size()==expectedGenomes && contigs.size()==expectedContigs, "Complete panel denominator changed");
	}
	private void legacy(String path, int expected, int expectedFull){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');int count=0, full=0;
		try{
			RrnaResourceIO.header(in, "locus_id\tgenome\tseqid\tstrand\tsource_start1\tsource_stop1\tcm_start1\tcm_stop1\tcm_state", path);
			for(byte[] row; (row=in.nextLine())!=null; ){
				p.set(row);require(p.terms()==9, "Accepted RF00174 panel ledger schema changed");
				final boolean supported=p.termEquals("CM_FULL_SPAN_SUPPORTED", 8);
				final int start=p.parseInt(6), stop=p.parseInt(7);
				require(supported?start>0 && stop>=start:start==0 && stop==0, "Legacy support and core coordinates disagree");
				add(p.parseString(0), p.parseString(1), p.parseString(2), p.parseString(3), 0,
					p.parseInt(4), p.parseInt(5), start, stop, supported?FAMILIES[0]:UNRESOLVED);
				count++;if(supported){full++;}
			}
		}finally{require(!in.close(), "Accepted RF00174 ledger read failed");}
		require(count==expected && full==expectedFull, "RF00174 truth must preserve its accepted denominators");
	}
	private void overlay(String path, int family, int expected){
		assert(family==1 || family==2):"Only sibling overlays use this schema; RF00174 retains its accepted panel ledger";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		final ObjectSet<String> seen=new ObjectSet<String>(String.class);int count=0;
		try{
			RrnaResourceIO.header(in, NcrnaClanOverlay.HEADER, path);
			for(byte[] row; (row=in.nextLine())!=null; ){
				p.set(row);require(p.terms()==33 && p.termEquals(FAMILIES[family], 21), "Sibling overlay family/schema changed");
				final String id=p.parseString(0);require(seen.add(id), "Duplicate overlay source ID: "+id);count++;
				final Genome g=genomes.get(p.parseString(1));if(g==null){continue;}
				require(p.termEquals(g.assembly, 2), "Same genome key must retain exact assembly identity: "+id);
				final String interval=p.parseString(8);final int dash=interval.indexOf('-');
				require(dash>0 && interval.indexOf(',')<0, "Reviewed overlay must retain contiguous source intervals: "+id);
				final int first=Integer.parseInt(interval.substring(0, dash)), last=Integer.parseInt(interval.substring(dash+1));
				require(last-first+1==p.parseInt(9), "Source interval length changed: "+id);
				final String member=p.parseString(31);final boolean admitted=p.termEquals("true", 23);
				int start=0, stop=0;
				if(admitted){
					start=p.parseInt(27);stop=p.parseInt(28);
					require(stop-start+1==p.parseInt(29) && (member.equals(FAMILIES[family]) && p.termEquals(member, 22)
						|| member.equals(CLAN_ONLY) && p.termEquals("CROSS_LABELED", 22)), "Reviewed biological labels/core length disagree: "+id);
				}else{require(p.termEquals("false", 23) && member.equals(UNRESOLVED) && p.termEquals(UNRESOLVED, 22), "Unadmitted rows must remain unresolved");}
				final Contig c=contigs.get(p.parseString(6));require(c!=null && c.length==p.parseInt(15), "Overlay contig length differs from actual panel FASTA: "+id);
				add(id, g.id, c.id, p.parseString(7), family, first, last, start, stop, member);
			}
		}finally{require(!in.close(), "Sibling overlay read failed");}
		require(count==expected, "Complete sibling overlay population changed: "+count+" expected="+expected);
	}
	private void add(String id, String genome, String seqid, String strand, int family, int first, int last, int start, int stop, String member){
		final Contig c=contigs.get(seqid);
		require(c!=null && c.genome.id.equals(genome) && (strand.equals("+") || strand.equals("-")), "Truth must bind its exact panel genome/contig/strand: "+id);
		require(first>0 && last>=first && last<=c.length && (start==0 && stop==0 || start>0 && stop>=start && stop<=c.length), "Truth coordinates leave their contig: "+id);
		require(ids.add(FAMILIES[family]+'/'+id), "Repeated source locus ID: "+id);
		require(physical.add(seqid+'\t'+strand+'\t'+first+'\t'+last), "Exact cross-source duplicate requires adjudication before union grading: "+id);
		if(start>0){require(cores.add(seqid+'\t'+strand+'\t'+start+'\t'+stop),
			"Distinct source intervals resolve to one exact reviewed core; adjudicate before FULL_CM union grading: "+seqid+':'+start+'-'+stop+" new="+id);}
		loci.add(new Locus(id, c, strand.charAt(0), family, first, last, start, stop, member));
	}

	/*--------------------------------------------------------------*/
	/*----------------            Output            ----------------*/
	/*--------------------------------------------------------------*/

	private void write(String out){
		RrnaResourceIO.fresh(out, ".loci.tsv", ".contigs.tsv", ".summary.tsv");
		final ByteStreamWriter rows=RrnaResourceIO.writer(out+".loci.tsv"), contigRows=RrnaResourceIO.writer(out+".contigs.tsv"), summary=RrnaResourceIO.writer(out+".summary.tsv");
		final ByteBuilder b=new ByteBuilder();final int[][] counts=new int[3][4];
		try{
			rows.println(HEADER);contigRows.println(CONTIG_HEADER);summary.println("source_family\thistorical\town_full\tclan_only_full\tunresolved");
			for(Genome g:orderedGenomes){for(Contig c:g.contigs){
				contigRows.print(b.clear().append(g.id).tab().append(g.assembly).tab().append(c.id).tab().append(c.length).nl());
			}}
			for(Locus l:loci){
				final int kind=l.member.equals(UNRESOLVED)?3:l.member.equals(CLAN_ONLY)?2:1;counts[l.family][0]++;counts[l.family][kind]++;
				rows.print(b.clear().append(l.id).tab().append(l.contig.genome.id).tab().append(l.contig.id).tab().append(FAMILIES[l.family]).tab().append(l.member)
					.tab().append(l.first).tab().append(l.last).tab().append(l.start).tab().append(l.stop).tab().append(l.strand).nl());
			}
			for(int i=0; i<3; i++){
				require(counts[i][0]==counts[i][1]+counts[i][2]+counts[i][3], "Each historical locus must retain exactly one evidence disposition");
				summary.print(b.clear().append(FAMILIES[i]).tab().append(counts[i][0]).tab().append(counts[i][1]).tab().append(counts[i][2]).tab().append(counts[i][3]).nl());
			}
		}finally{RrnaResourceIO.finish(rows, contigRows, summary);}
		System.err.println("COBALAMIN_PANEL_PREP_PASS genomes="+genomes.size()+" contigs="+contigs.size()+" loci="+loci.size());
	}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final class Genome {
		Genome(String id_, String assembly_, String fasta_){id=id_;assembly=assembly_;fasta=fasta_;}
		final String id, assembly, fasta;long bases;final ArrayList<Contig> contigs=new ArrayList<Contig>();
	}
	private static final class Contig {
		Contig(Genome g, String id_, int length_){genome=g;id=id_;length=length_;}
		final Genome genome;final String id;final int length;
	}
	private static final class Locus {
		Locus(String id_, Contig c, char strand_, int family_, int first_, int last_, int start_, int stop_, String member_){
			id=id_;contig=c;strand=strand_;family=family_;first=first_;last=last_;start=start_;stop=stop_;member=member_;
		}
		final String id, member;final Contig contig;final char strand;final int family, first, last, start, stop;
	}
	private final ObjectMap<String,Genome> genomes=new ObjectMap<String,Genome>(String.class, Genome.class);
	private final ObjectMap<String,Contig> contigs=new ObjectMap<String,Contig>(String.class, Contig.class);
	private final ObjectSet<String> ids=new ObjectSet<String>(String.class), physical=new ObjectSet<String>(String.class);
	private final ObjectSet<String> cores=new ObjectSet<String>(String.class);
	private final ArrayList<Genome> orderedGenomes=new ArrayList<Genome>();
	private final ArrayList<Locus> loci=new ArrayList<Locus>();
	public static final String[] FAMILIES={"RF00174", "RF01482", "RF01689"};
	public static final String CLAN_ONLY="B12clade", UNRESOLVED="UNRESOLVED";
	public static final String HEADER="locus_id\tgenome\tseqid\tsource_family\tmember_family\tsource_start1\tsource_stop1\tcore_start1\tcore_stop1\tstrand";
	public static final String CONTIG_HEADER="genome\tassembly\tseqid\tlength";
}
