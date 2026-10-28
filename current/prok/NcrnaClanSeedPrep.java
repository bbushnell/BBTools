package prok;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectMap;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Joins complete source corpora and admitted corrected cores before one clan split.
 * Raw annotation rows remain unchanged; corrected seed geometry has a separate schema.
 * @author Brian, Raiden
 */
public final class NcrnaClanSeedPrep {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args, NcrnaClanSeedPrep.class, false);
		String in=null, manifest=null, out=null;int expected=-1, expectedGenomes=-1, expectedSource=-1;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0, "Use flag=value");
			final String k=arg.substring(0, eq).toLowerCase(Locale.ROOT), v=arg.substring(eq+1);
			if(k.equals("in")){in=v;}else if(k.equals("cohorts")){manifest=v;}else if(k.equals("outprefix")){out=v;}
			else if(k.equals("expected")){expected=Integer.parseInt(v);}else if(k.equals("expectedgenomes")){expectedGenomes=Integer.parseInt(v);}
			else if(k.equals("expectedsource")){expectedSource=Integer.parseInt(v);}else{throw new IllegalArgumentException("Unknown option: "+k);}
		}
		require(in!=null && manifest!=null && out!=null && expected>0 && expectedGenomes>0 && expectedSource>=expected,
			"Bind the complete admitted FASTA, three full source cohorts and explicit population denominators");
		Shared.setThreads(1);Shared.TRIM_READ_DESCRIPTION=false;
		try{final NcrnaClanSeedPrep x=new NcrnaClanSeedPrep();x.load(manifest, expectedGenomes, expectedSource);x.select(in, expected);x.write(out);}
		finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	private void load(String path, int expectedGenomes, int expectedSource){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{
			RrnaResourceIO.header(in, "family\tloci\tgenomes\thashes\texpected", path);
			for(byte[] row; (row=in.nextLine())!=null; ){
				p.set(row);final int i=cohorts.size();require(p.terms()==5 && i<3 && p.termEquals(FAMILIES[i], 0), "Declare all three source families in canonical order");
				final Cohort c=new Cohort(i, p);require(c.expected>0 && c.loci.startsWith("/") && c.genomes.startsWith("/") && c.hashes.startsWith("/"),
					"Source cohorts require positive denominators and absolute original paths");cohorts.add(c);
			}
		}finally{require(!in.close(), "Cohort manifest read failed");}
		require(cohorts.size()==3, "The clan must include all three declared source families");
		for(Cohort c:cohorts){genomes(c, expectedGenomes);hashes(c.hashes);loci(c);}
		require(raw.size()==expectedSource, "Combined raw annotation denominator changed: "+raw.size()+" expected="+expectedSource);
		for(Genome g:genomes){
			require(g.mask==7, "Every source cohort must contain the same complete genome roster: "+g.name);
			int total=0;for(int i=0; i<3; i++){require(g.declared[i]==g.observed[i], "Raw locus/genome counts disagree for "+g.name+" family="+FAMILIES[i]);total=Math.addExact(total, g.declared[i]);}
			require(total==0 || hashes.contains(g.fasta), "Every positive genome in the combined clan needs a bound FASTA hash: "+g.name);
			g.total=total;
		}
	}
	private void genomes(Cohort c, int expected){
		final ByteFile in=RrnaResourceIO.open(c.genomes);final LineParser1 p=new LineParser1('\t'), old=new LineParser1('\t');int count=0;
		try{
			RrnaResourceIO.header(in, NcrnaCorpusSample.GENOME_HEADER, c.genomes);
			for(byte[] row; (row=in.nextLine())!=null; ){
				p.set(row);require(p.terms()==11 && p.parseInt(8)>=0, "Complete genome rows need the canonical schema and nonnegative locus counts");
				final String name=p.parseString(1);Genome g=byGenome.get(name);
				if(c.index==0){require(g==null, "Repeated source genome: "+name);g=new Genome(name, p.parseString(7), row);byGenome.put(name, g);genomes.add(g);}
				else{
					require(g!=null, "Sibling cohort introduced a different genome universe: "+name);old.set(g.row);
					for(int i=0; i<8; i++){require(same(p, old, i, i), "Shared genome identity/path metadata differs: "+name+" column="+i);}
				}
				final int bit=1<<c.index;require((g.mask&bit)==0, "Repeated genome within one family cohort: "+name);
				g.mask|=bit;g.declared[c.index]=p.parseInt(8);count++;
				final long bases=p.parseLong(9);final int contigs=p.parseInt(10);
				// NcrnaAnnotatedCorpus.extractGenome returns before reading a negative
				// genome, so its bases/contigs are0, not a conflicting size measurement.
				require(g.declared[c.index]==0?bases==0 && contigs==0:bases>0 && contigs>0,
					"Genome size observations must reflect whether this family read the FASTA: "+name);
				if(bases>0){require(g.bases==0 || g.bases==bases && g.contigs==contigs, "Positive-family FASTA size observations disagree: "+name);g.bases=bases;g.contigs=contigs;}
			}
		}finally{require(!in.close(), "Genome roster read failed");}
		require(count==expected && genomes.size()==expected, "Every cohort must preserve the complete declared genome population");
	}
	private void hashes(String path){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{
			for(byte[] row; (row=in.nextLine())!=null; ){
				p.set(row);require(p.terms()==2, "Source provenance must contain sha80 and absolute path");
				final String digest=p.parseString(0), file=p.parseString(1), prior=hashes.get(file);
				require(digest.matches("[0-9a-f]{20}") && file.startsWith("/") && (prior==null || prior.equals(digest)),
					"Shared paths must have consistent sha80 provenance across the three frozen corpora: "+file);
				if(prior==null){hashes.put(file, digest);hashPaths.add(file);}
			}
		}finally{require(!in.close(), "Source hash ledger read failed");}
	}
	private void loci(Cohort c){
		final ByteFile in=RrnaResourceIO.open(c.loci);final LineParser1 p=new LineParser1('\t'), g=new LineParser1('\t');int count=0;
		try{
			RrnaResourceIO.header(in, NcrnaCmLocusAudit.HEADER, c.loci);
			for(byte[] row; (row=in.nextLine())!=null; ){
				p.set(row);require(p.terms()==21, "Raw annotation schema changed");final String id=p.parseString(0);
				final Genome genome=byGenome.get(p.parseString(1));
				require(id.startsWith(FAMILIES[c.index]+"_") && !byLocus.contains(id) && genome!=null && p.parseInt(9)>0 && p.parseInt(15)>0,
					"Every raw source locus needs a unique family identity, known genome and positive lengths: "+id);
				g.set(genome.row);for(int i=2; i<=5; i++){require(same(p, g, i, i), "Source locus and genome disagree on assembly/taxonomy: "+id);}
				final Locus l=new Locus(id, row, genome, c.index);byLocus.put(id, l);raw.add(l);genome.observed[c.index]++;count++;
			}
		}finally{require(!in.close(), "Raw annotation read failed");}
		require(count==c.expected, "Full source family population changed: "+FAMILIES[c.index]);
	}
	private void select(String fasta, int expected){
		final LineParser1 p=new LineParser1('\t');
		RrnaResourceIO.stream(fasta, r->{
			final String id=RrnaResourceIO.first(r.id);final Locus l=byLocus.get(id);
			require(l!=null && !l.selected && r.bases!=null && r.length()>=17, "Every admitted member must map once to its original source: "+id);
			p.set(l.row);require(RrnaResourceIO.field(r.id, "source").equals(l.genome.name), "Admitted header source disagrees with original annotation: "+id);
			final String seqid=RrnaResourceIO.optionalField(r.id, "contig"), strand=RrnaResourceIO.optionalField(r.id, "strand"), sourceFamily=RrnaResourceIO.optionalField(r.id, "source_family");
			require((seqid==null || p.termEquals(seqid, 6)) && (strand==null || p.termEquals(strand, 7)) && (sourceFamily==null || sourceFamily.equals(FAMILIES[l.family])),
				"Corrected core must preserve original contig, strand and source-family identity: "+id);
			final String a=RrnaResourceIO.optionalField(r.id, "start"), z=RrnaResourceIO.optionalField(r.id, "stop");
			require((a==null)==(z==null), "Corrected bounds must supply both start and stop: "+id);
			final int start, stop;
			if(a==null){
				final String interval=p.parseString(8);final int dash=interval.indexOf('-');
				require(dash>0 && interval.indexOf(',')<0 && interval.indexOf('-', dash+1)<0, "Uncorrected admitted inputs must retain one contiguous original interval: "+id);
				start=Integer.parseInt(interval.substring(0, dash));stop=Integer.parseInt(interval.substring(dash+1));
			}else{start=Integer.parseInt(a);stop=Integer.parseInt(z);}
			require(start>0 && stop<=p.parseInt(15) && (long)stop-start+1==r.length(), "Seed alias coordinates must describe the actual admitted core length within its contig: "+id);
			for(byte base:r.bases){require(RrnaResourceIO.base(base)>=0, "Admitted seed cores must contain unambiguous bases: "+id);}
			final String biological=RrnaResourceIO.field(r.id, "member_family");final int group=index(RrnaResourceIO.field(r.id, "build_group"));
			final int member=biological.equals("B12clade")?3:index(biological);
			require(group>=0 && (member==l.family && group==l.family || member==3), "Biological family and operational group must retain the reviewed own-or-clan policy: "+id);
			l.selected=true;selected.add(new Seed(l, start, stop, member, group));admittedSource[l.family]++;admittedBiological[member]++;buildGroups[group]++;
		});
		require(selected.size()==expected, "The complete admitted population changed: "+selected.size()+" expected="+expected);
	}
	private void write(String out){
		RrnaResourceIO.fresh(out, ".raw.loci.tsv.gz", ".genomes.tsv", ".hashes.sha80", ".seed.loci.tsv", ".eligible.ids.txt", ".summary.tsv");
		final ByteStreamWriter loci=RrnaResourceIO.writer(out+".raw.loci.tsv.gz"), roster=RrnaResourceIO.writer(out+".genomes.tsv"), hash=RrnaResourceIO.writer(out+".hashes.sha80");
		final ByteStreamWriter seeds=RrnaResourceIO.writer(out+".seed.loci.tsv"), ids=RrnaResourceIO.writer(out+".eligible.ids.txt"), summary=RrnaResourceIO.writer(out+".summary.tsv");
		final LineParser1 p=new LineParser1('\t');final ByteBuilder b=new ByteBuilder();
		try{
			loci.println(NcrnaCmLocusAudit.HEADER);for(Locus l:raw){loci.println(l.row);}
			roster.println(NcrnaCorpusSample.GENOME_HEADER);
			for(Genome g:genomes){
				p.set(g.row);b.clear();for(int i=0; i<11; i++){if(i>0){b.tab();}if(i==8){b.append(g.total);}else if(i==9){b.append(g.bases);}else if(i==10){b.append(g.contigs);}else{p.appendTerm(b, i);}}roster.print(b.nl());
			}
			Collections.sort(hashPaths);for(String file:hashPaths){hash.print(b.clear().append(hashes.get(file)).tab().append(file).nl());}
			seeds.println(SEED_HEADER);
			for(Seed s:selected){
				ids.println(s.locus.id);p.set(s.locus.row);b.clear();
				for(int i=0; i<8; i++){if(i>0){b.tab();}p.appendTerm(b, i);}
				seeds.print(b.tab().append(s.start).append('-').append(s.stop).tab().append(FAMILIES[s.locus.family]).tab()
					.append(s.member==3?"B12clade":FAMILIES[s.member]).tab().append(FAMILIES[s.group]).nl());
			}
			summary.println("view\tlabel\tcount");
			for(int i=0; i<3; i++){
				summary.print(b.clear().append("RAW_SOURCE\t").append(FAMILIES[i]).tab().append(cohorts.get(i).expected).nl());
				summary.print(b.clear().append("ADMITTED_SOURCE\t").append(FAMILIES[i]).tab().append(admittedSource[i]).nl());
				summary.print(b.clear().append("BUILD_GROUP\t").append(FAMILIES[i]).tab().append(buildGroups[i]).nl());
			}
			for(int i=0; i<4; i++){summary.print(b.clear().append("ADMITTED_BIOLOGICAL\t").append(i==3?"B12clade":FAMILIES[i]).tab().append(admittedBiological[i]).nl());}
			summary.print(b.clear().append("TOTAL\tRAW_LOCI\t").append(raw.size()).append("\nTOTAL\tADMITTED_LOCI\t").append(selected.size())
				.append("\nTOTAL\tGENOMES\t").append(genomes.size()).append("\nTOTAL\tHASHED_PATHS\t").append(hashPaths.size()).nl());
		}finally{RrnaResourceIO.finish(loci, roster, hash, seeds, ids, summary);}
		System.err.println("NCRNA_CLAN_SEED_PREP_PASS source="+raw.size()+" admitted="+selected.size()+" genomes="+genomes.size()+"; raw annotations preserved, corrected seed geometry separate");
	}
	private static boolean same(LineParser1 a, LineParser1 b, int x, int y){
		a.setBounds(x);b.setBounds(y);final int n=a.b()-a.a();if(n!=b.b()-b.a()){return false;}
		final byte[] aa=(byte[])a.line(), bb=(byte[])b.line();for(int i=0; i<n; i++){if(aa[a.a()+i]!=bb[b.a()+i]){return false;}}return true;
	}
	private static int index(String family){for(int i=0; i<3; i++){if(family.equals(FAMILIES[i])){return i;}}return -1;}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final class Cohort {
		Cohort(int i, LineParser1 p){index=i;loci=p.parseString(1);genomes=p.parseString(2);hashes=p.parseString(3);expected=p.parseInt(4);}
		final int index, expected;final String loci, genomes, hashes;
	}
	private static final class Genome {
		Genome(String n, String f, byte[] r){name=n;fasta=f;row=r;}
		final String name, fasta;final byte[] row;final int[] declared=new int[3], observed=new int[3];int mask, total, contigs;long bases;
	}
	private static final class Locus {
		Locus(String i, byte[] r, Genome g, int f){id=i;row=r;genome=g;family=f;}
		final String id;final byte[] row;final Genome genome;final int family;boolean selected;
	}
	private static final class Seed {
		Seed(Locus l, int a, int z, int m, int g){locus=l;start=a;stop=z;member=m;group=g;}
		final Locus locus;final int start, stop, member, group;
	}
	private final ArrayList<Cohort> cohorts=new ArrayList<Cohort>();private final ArrayList<Genome> genomes=new ArrayList<Genome>();
	private final ArrayList<Locus> raw=new ArrayList<Locus>();private final ArrayList<Seed> selected=new ArrayList<Seed>();private final ArrayList<String> hashPaths=new ArrayList<String>();
	private final ObjectMap<String,Genome> byGenome=new ObjectMap<String,Genome>(String.class, Genome.class);
	private final ObjectMap<String,Locus> byLocus=new ObjectMap<String,Locus>(String.class, Locus.class);
	private final ObjectMap<String,String> hashes=new ObjectMap<String,String>(String.class, String.class);
	private final int[] admittedSource=new int[3], admittedBiological=new int[4], buildGroups=new int[3];
	private static final String[] FAMILIES={"RF00174", "RF01482", "RF01689"};
	public static final String SEED_HEADER="locus_id\tgenome\tassembly\ttaxid\tdomain\tphylum\tseqid\tstrand\tintervals1\tsource_family\tmember_family\tbuild_group";
}
