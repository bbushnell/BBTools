package prok;

import java.util.ArrayList;
import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Projects reviewed B12 evidence onto one exact genome panel without changing biological labels.
 * Optional maincorpus, mainrows, rf14roster and rf16roster inputs expand the
 * original panel to every admitted rare-subtype/clan-only member's genome.
 * Additional main-subtype annotations are UNREVIEWED, never inferred full cores.
 * In this mode, genomes/contigs still bind the original roster; inputs.tsv records
 * the expanded roster and summary.tsv adds an unreviewed column.
 * @author Brian Bushnell, Raiden
 */
public final class CobalaminPanelPrep {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args, CobalaminPanelPrep.class, false);
		String roster=null, legacy=null, rf14=null, rf16=null, out=null;
		String mainCorpus=null, roster14=null, roster16=null;
		int genomes=-1, contigs=-1, oldRows=-1, oldFull=-1, rows14=-1, rows16=-1;
		int mainRows=-1;
		for(String arg:pp.args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT);
			final String b=split.length>1 ? split[1] : null;
			require(!a.isEmpty() && b!=null, "Use flag=value: "+arg);
			if(a.equals("roster")){roster=b;}
			else if(a.equals("legacy")){legacy=b;}
			else if(a.equals("rf14")){rf14=b;}
			else if(a.equals("rf16")){rf16=b;}
			else if(a.equals("outprefix")){out=b;}
			else if(a.equals("genomes")){genomes=Integer.parseInt(b);}
			else if(a.equals("contigs")){contigs=Integer.parseInt(b);}
			else if(a.equals("legacyrows")){oldRows=Integer.parseInt(b);}
			else if(a.equals("legacyfull")){oldFull=Integer.parseInt(b);}
			else if(a.equals("rf14rows")){rows14=Integer.parseInt(b);}
			else if(a.equals("rf16rows")){rows16=Integer.parseInt(b);}
			else if(a.equals("maincorpus")){mainCorpus=b;}
			else if(a.equals("mainrows")){mainRows=Integer.parseInt(b);}
			else if(a.equals("rf14roster")){roster14=b;}
			else if(a.equals("rf16roster")){roster16=b;}
			else{throw new IllegalArgumentException("Unknown option: "+a);}
		}
		require(roster!=null && legacy!=null && rf14!=null && rf16!=null && out!=null && genomes>0 && contigs>0
			&& oldRows>=0 && oldFull>=0 && rows14>=0 && rows16>=0, "Bind the exact panel and complete reviewed evidence populations");
		final boolean expanded=mainCorpus!=null || roster14!=null || roster16!=null || mainRows>=0;
		require(!expanded || mainCorpus!=null && roster14!=null && roster16!=null && mainRows>=oldRows,
			"Expansion requires the complete main corpus/count and both exact rare-subtype rosters");
		Shared.setThreads(1);Shared.TRIM_READ_DESCRIPTION=false;
		try{
			final CobalaminPanelPrep x=new CobalaminPanelPrep();
			x.expanded=expanded;
			x.roster(roster, genomes, contigs);
			if(expanded){
				x.expandRoster(roster14, rf14, 1);
				x.expandRoster(roster16, rf16, 2);
			}
			x.legacy(legacy, oldRows, oldFull);
			if(expanded){x.mainCorpus(mainCorpus, mainRows);}
			x.overlay(rf14, 1, rows14);
			x.overlay(rf16, 2, rows16);
			x.write(out);
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
				addGenome(p, row);
			}
		}finally{require(!in.close(), "Panel roster read failed");}
		require(genomes.size()==expectedGenomes && contigs.size()==expectedContigs, "Complete panel denominator changed");
	}

	/** Reads a selected whole genome once, preserving its metadata and exact contig ownership. */
	private void addGenome(LineParser1 p, byte[] row){
		final Genome g=new Genome(p.parseString(1), p.parseString(2), p.parseString(7), row);
		require(!genomes.contains(g.id), "Duplicate panel genome: "+g.id);
		genomes.put(g.id, g);
		orderedGenomes.add(g);
		final long expectedBases=p.parseLong(9);
		final int expectedN=p.parseInt(10);
		RrnaResourceIO.stream(g.fasta, r->{
			final String name=RrnaResourceIO.first(r.id);
			require(r.bases!=null && r.length()>0 && !contigs.contains(name),
				"Combined grading requires unique nonempty contigs: "+name);
			final Contig c=new Contig(g, name, r.length());
			contigs.put(name, c);
			g.contigs.add(c);
			g.bases+=r.length();
		});
		// Original extraction may leave a target-free genome's expected counts at zero.
		require(expectedBases==0 && expectedN==0 || g.bases==expectedBases && g.contigs.size()==expectedN,
			"Measured source genome bases/contigs changed: "+g.id);
	}

	/** Adds genomes named by admitted overlay members, including clan-only evidence. */
	private void expandRoster(String roster, String overlay, int family){
		assert(family==1 || family==2):"Expansion is sourced only from the two reviewed sibling overlays";
		final ObjectMap<String,String> selected=new ObjectMap<String,String>(String.class, String.class);
		final LineParser1 p=new LineParser1('\t');
		final ByteFile evidence=RrnaResourceIO.open(overlay);
		try{
			RrnaResourceIO.header(evidence, NcrnaClanOverlay.HEADER, overlay);
			for(byte[] row; (row=evidence.nextLine())!=null; ){
				p.set(row);
				require(p.terms()==33 && p.termEquals(FAMILIES[family], 21), "Expansion overlay family/schema changed");
				if(!p.termEquals("true", 23)){continue;}
				final String id=p.parseString(1), assembly=p.parseString(2);
				final String prior=selected.get(id);
				require(prior==null || prior.equals(assembly), "Overlay genome key has conflicting assemblies: "+id);
				selected.put(id, assembly);
			}
		}finally{require(!evidence.close(), "Expansion overlay read failed");}
		final ObjectSet<String> found=new ObjectSet<String>(String.class);
		final ByteFile in=RrnaResourceIO.open(roster);
		try{
			RrnaResourceIO.header(in, NcrnaCalibrationPrep.ROSTER_HEADER, roster);
			for(byte[] row; (row=in.nextLine())!=null; ){
				p.set(row);
				require(p.terms()==11, "Expansion roster schema changed");
				final String id=p.parseString(1), assembly=selected.get(id);
				if(assembly==null){continue;}
				require(found.add(id) && p.termEquals(assembly, 2), "Admitted genome needs one exact roster assembly: "+id);
				final Genome existing=genomes.get(id);
				if(existing==null){addGenome(p, row);}
				else{
					require(existing.assembly.equals(assembly) && p.termEquals(existing.fasta, 7),
						"Shared genome must retain its exact assembly and FASTA path: "+id);
				}
			}
		}finally{require(!in.close(), "Expansion roster read failed");}
		require(found.size()==selected.size(), "Every admitted sibling member needs its genome in the supplied roster");
	}

	/** Preserves reviewed main loci verbatim and masks other selected SOURCE annotations. */
	private void mainCorpus(String path, int expected){
		final ObjectMap<String,Locus> reviewed=new ObjectMap<String,Locus>(String.class, Locus.class);
		for(Locus l:loci){
			assert(l.family==0):"mainCorpus must follow legacy loading and precede sibling overlays";
			reviewed.put(l.id, l);
		}
		final ObjectSet<String> seen=new ObjectSet<String>(String.class);
		final ByteFile in=RrnaResourceIO.open(path);
		final LineParser1 p=new LineParser1('\t');
		int count=0, retained=0;
		try{
			RrnaResourceIO.header(in, NcrnaCmLocusAudit.HEADER, path);
			for(byte[] row; (row=in.nextLine())!=null; ){
				p.set(row);
				require(p.terms()==21, "Main SOURCE ledger must retain all original source/context columns");
				count++;
				final Genome g=genomes.get(p.parseString(1));
				if(g==null){continue;}
				final String id=p.parseString(0), seqid=p.parseString(6), strand=p.parseString(7);
				require(seen.add(id) && p.termEquals(g.assembly, 2), "Main SOURCE identity or assembly changed: "+id);
				final String interval=p.parseString(8);
				final int dash=interval.indexOf('-');
				require(dash>0 && interval.indexOf(',')<0, "Main SOURCE interval needs contiguous reviewed geometry: "+id);
				final int first=Integer.parseInt(interval.substring(0, dash));
				final int last=Integer.parseInt(interval.substring(dash+1));
				final Contig c=contigs.get(seqid);
				require(c!=null && c.genome==g && c.length==p.parseInt(15) && last-first+1==p.parseInt(9),
					"Main SOURCE coordinates and ownership must match the actual FASTA: "+id);
				final Locus prior=reviewed.get(id);
				if(prior==null){add(id, g.id, seqid, strand, 0, first, last, 0, 0, UNREVIEWED);}
				else{
					require(prior.contig==c && prior.first==first && prior.last==last && strand.length()==1 && strand.charAt(0)==prior.strand,
						"Legacy SOURCE coordinates must agree before preserving the existing review: "+id);
					retained++;
				}
			}
		}finally{require(!in.close(), "Main SOURCE ledger read failed");}
		require(count==expected && retained==reviewed.size(), "Main corpus must retain every legacy locus and its complete input population");
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
		c.genome.locusCount++;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Output            ----------------*/
	/*--------------------------------------------------------------*/

	private void write(String out){
		RrnaResourceIO.fresh(out, ".loci.tsv", ".contigs.tsv", ".summary.tsv");
		if(expanded){RrnaResourceIO.fresh(out, ".inputs.tsv");}
		final ByteStreamWriter rows=RrnaResourceIO.writer(out+".loci.tsv"), contigRows=RrnaResourceIO.writer(out+".contigs.tsv"), summary=RrnaResourceIO.writer(out+".summary.tsv");
		final ByteBuilder b=new ByteBuilder();final int[][] counts=new int[3][5];
		try{
			rows.println(HEADER);
			contigRows.println(CONTIG_HEADER);
			summary.println("source_family\thistorical\town_full\tclan_only_full\tunresolved"+(expanded ? "\tunreviewed" : ""));
			for(Genome g:orderedGenomes){for(Contig c:g.contigs){
				contigRows.print(b.clear().append(g.id).tab().append(g.assembly).tab().append(c.id).tab().append(c.length).nl());
			}}
			for(Locus l:loci){
				final int kind=l.member.equals(UNREVIEWED) ? 4 : l.member.equals(UNRESOLVED) ? 3 : l.member.equals(CLAN_ONLY) ? 2 : 1;
				counts[l.family][0]++;
				counts[l.family][kind]++;
				rows.print(b.clear().append(l.id).tab().append(l.contig.genome.id).tab().append(l.contig.id).tab().append(FAMILIES[l.family]).tab().append(l.member)
					.tab().append(l.first).tab().append(l.last).tab().append(l.start).tab().append(l.stop).tab().append(l.strand).nl());
			}
			for(int i=0; i<3; i++){
				require(counts[i][0]==counts[i][1]+counts[i][2]+counts[i][3]+counts[i][4],
					"Each historical locus must retain exactly one evidence disposition");
				b.clear().append(FAMILIES[i]).tab().append(counts[i][0]).tab().append(counts[i][1]).tab().append(counts[i][2]).tab().append(counts[i][3]);
				if(expanded){b.tab().append(counts[i][4]);}
				summary.print(b.nl());
			}
		}finally{RrnaResourceIO.finish(rows, contigRows, summary);}
		if(expanded){writeRoster(out+".inputs.tsv");}
		System.err.println("COBALAMIN_PANEL_PREP_PASS genomes="+genomes.size()+" contigs="+contigs.size()+" loci="+loci.size());
	}

	/** Writes one caller input per selected genome with measured whole-FASTA sizes. */
	private void writeRoster(String path){
		assert(expanded):"Only expanded panels add a roster; the legacy output contract remains unchanged";
		final ByteStreamWriter out=RrnaResourceIO.writer(path);
		final ByteBuilder b=new ByteBuilder();
		final LineParser1 p=new LineParser1('\t');
		try{
			out.println(NcrnaCalibrationPrep.ROSTER_HEADER);
			for(int i=0; i<orderedGenomes.size(); i++){
				final Genome g=orderedGenomes.get(i);
				p.set(g.roster);
				b.clear().append(i);
				for(int term=1; term<=7; term++){p.appendTerm(b.tab(), term);}
				out.print(b.tab().append(g.locusCount).tab().append(g.bases).tab().append(g.contigs.size()).nl());
			}
		}finally{RrnaResourceIO.finish(out);}
	}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final class Genome {
		Genome(String id_, String assembly_, String fasta_, byte[] roster_){
			id=id_;
			assembly=assembly_;
			fasta=fasta_;
			roster=roster_.clone();
		}
		final String id, assembly, fasta;long bases;final ArrayList<Contig> contigs=new ArrayList<Contig>();
		final byte[] roster;
		int locusCount;
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
	/** Additional corpus/roster inputs enable expansion; omitted inputs preserve legacy behavior. */
	private boolean expanded=false;
	public static final String[] FAMILIES={"RF00174", "RF01482", "RF01689"};
	public static final String CLAN_ONLY="B12clade", UNRESOLVED="UNRESOLVED", UNREVIEWED="UNREVIEWED";
	public static final String HEADER="locus_id\tgenome\tseqid\tsource_family\tmember_family\tsource_start1\tsource_stop1\tcore_start1\tcore_stop1\tstrand";
	public static final String CONTIG_HEADER="genome\tassembly\tseqid\tlength";
}
