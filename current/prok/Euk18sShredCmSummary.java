package prok;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Descriptive primary-CM evidence per shred, never automatic truth adjudication.
 * Keeps no-hit targets in denominators and preserves ties/partial/truncation flags.
 * Raw tblout remains the source for all hits, including non-primary loci.
 * fullModelHits/primaryFullModel require inclusion, model span 1..CLEN AND trunc=no;
 * they are not pure model-span coverage. Model coordinates and truncation remain separate.
 * Equal scores use lexical retained-field order for reproducibility, not biological preference.
 * @author Raiden
 */
public final class Euk18sShredCmSummary {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args, Euk18sShredCmSummary.class, false);String regions=null, cm=null, table=null, out=null, accession=null;int flank=-1;
		for(String arg:pp.args){final int e=arg.indexOf('=');require(e>0, "Expected flag=value");final String k=arg.substring(0, e).toLowerCase(Locale.ROOT), v=arg.substring(e+1);
			if(k.equals("regions")){regions=v;}else if(k.equals("cm")){cm=v;}else if(k.equals("tblout")){table=v;}else if(k.equals("outprefix")){out=v;}
			else if(k.equals("accession")){accession=v;}else if(k.equals("flank")){flank=Integer.parseInt(v);}else{throw new IllegalArgumentException(k);}}
		require(regions!=null && cm!=null && table!=null && out!=null && accession!=null && flank>=0, "Declare regions/cm/tblout/outprefix/accession/flank");
		Shared.setThreads(1);new Euk18sShredCmSummary(accession, flank).run(regions, cm, table, out);Shared.closeStream(pp.outstream);
	}
	public Euk18sShredCmSummary(String acc, int f){require(acc.matches("RF[0-9]{5}") && f>=0, "One explicit Rfam model and flank required");accession=acc;flank=f;}
	public void run(String regions, String cm, String table, String prefix){
		RrnaResourceIO.fresh(prefix, ".primary.tsv", ".summary.tsv");loadModel(cm);loadRegions(regions);loadHits(table);write(prefix);
		System.err.println("SHRED_CM_SUMMARY_PASS model="+accession+" clen="+clen+" targets="+ordered.size()+" rawHits="+rawHits+" newAlignments=0");
	}
	void loadModel(String path){
		final ByteFile in=RrnaResourceIO.open(path);String actual=null;
		try{final byte[] header=in.nextLine();require(header!=null && text(header).startsWith("INFERNAL1/"), "Expected one Infernal CM file");
			for(byte[] row; (row=in.nextLine())!=null;){
				if(starts(row, "ACC ")){actual=text(row).substring(4).trim();}
				else if(starts(row, "CLEN ")){clen=Integer.parseInt(text(row).substring(5).trim());}
				if(actual!=null && clen>0){break;}
			}
		}finally{require(!in.close(), "CM header read failed");}require(accession.equals(actual) && clen>0, "Actual CM accession/CLEN must match the declared model");
	}
	void loadRegions(String path){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{RrnaResourceIO.header(in, Euk18sShredSet.REGION_HEADER, path);for(byte[] row; (row=in.nextLine())!=null;){p.set(row);require(p.terms()==13, "Region schema changed");if(p.parseInt(1)!=flank){continue;}
			final int a=p.parseInt(6), z=p.parseInt(7);final long n=(long)z-a+1;require(a>0 && n>0 && n<=Integer.MAX_VALUE, "Shred length requires valid inclusive source coordinates");
			final Target t=new Target(p.parseString(0), p.parseString(2), p.parseString(3), p.parseString(4), (int)n);Euk18sShredSet.roleIndex(t.role);
			require(byId.put(t.id, t)==null, "Duplicate shred would inflate the scan denominator");ordered.add(t);
		}}finally{require(!in.close(), "Region read failed");}require(!ordered.isEmpty(), "No targets in the requested flank");
	}
	void loadHits(String path){
		final ByteFile in=RrnaResourceIO.open(path);final RrnaCmTruthBuilder.Fields f=new RrnaCmTruthBuilder.Fields();boolean program=false, version=false, done=false;
		try{for(byte[] row; (row=in.nextLine())!=null;){if(row.length==0){continue;}require(!done, "Content after Infernal completion marker");
			if(row[0]=='#'){final String s=text(row);if(s.startsWith("# Program:")){require(s.substring(10).trim().equals("cmsearch"), "Expected cmsearch output");program=true;}
				if(s.startsWith("# Version:")){require(s.substring(10).trim().startsWith("1.1.5 "), "Expected the pinned Infernal1.1.5 table format");version=true;}
				if(s.equals("# [ok]")){done=true;}continue;}
			f.set(row);require(f.eq(3, accession) && f.eq(4, "cm"), "Tblout model differs from declared accession; never relabel another model's hits");
			final Target t=byId.get(f.text(0));require(t!=null, "CM target absent from the bound shred manifest");
			final int mf=f.number(5), mt=f.number(6), from=f.number(7), to=f.number(8);final String strand=f.text(9), trunc=f.text(10);
			require(mf>=1 && mt>=mf && mt<=clen && from>=1 && to>=1 && from<=t.length && to<=t.length, "CM model/sequence coordinates exceed their actual bounds");
			require((strand.equals("+") && from<=to) || (strand.equals("-") && from>=to), "CM strand must agree with sequence coordinate order");
			require(trunc.equals("no") || trunc.equals("5'") || trunc.equals("3'") || trunc.equals("5'&3'"), "Unknown CM truncation flag");
			require(f.eq(16, "!") || f.eq(16, "?"), "Unknown CM inclusion marker");final boolean included=f.eq(16, "!");
			final float score=Parse.parseFloat(f.bytes, f.start[14], f.end[14]);
			// E-values use scientific notation; retain exponent-aware parsing explicitly.
			final double evalue=Double.parseDouble(f.text(15));require(Float.isFinite(score) && Double.isFinite(evalue) && evalue>=0, "CM score/E-value must be finite and E-value nonnegative");
			final boolean full=included && mf==1 && mt==clen && trunc.equals("no");t.hits++;rawHits++;if(included){t.included++;}if(full){t.full++;}
			final int a=Math.min(from, to), z=Math.max(from, to);final char sign=strand.charAt(0);
			if(t.best==null || score>t.best.score){t.best=new Hit(mf, mt, a, z, sign, trunc, included, full, score, evalue);t.ties=1;}
			else if(score==t.best.score){t.ties++;if(precedes(t.best, mf, mt, a, z, sign, trunc, included, score, evalue)){t.best=new Hit(mf, mt, a, z, sign, trunc, included, full, score, evalue);}}
		}}finally{require(!in.close(), "CM table read failed");}require(program && version && done, "Missing cmsearch/version/completion footer");
	}
	void write(String prefix){
		final ByteStreamWriter primary=RrnaResourceIO.writer(prefix+".primary.tsv"), summary=RrnaResourceIO.writer(prefix+".summary.tsv");final ByteBuilder b=new ByteBuilder();
		try{primary.println("shred\tclade\tassembly\trole\tbases\treportedHits\tincludedHits\tfullModelHits\tprimaryModel\tprimaryBits\tprimaryEvalue\tbestScoreTies\tmodelFrom\tmodelTo\tclen\tseqStart1\tseqStop1\tstrand\ttrunc\tincluded\tprimaryFullModel\tprimaryAtLeastHalfLength");
			for(Target t:ordered){b.clear().append(t.id).tab().append(t.clade).tab().append(t.assembly).tab().append(t.role).tab().append(t.length).tab().append(t.hits).tab().append(t.included).tab().append(t.full);
				final Hit h=t.best;if(h==null){for(int i=0; i<14; i++){b.append("\t.");}}
				else{b.tab().append(accession).tab().appendSlow(h.score).tab().appendSlow(h.evalue).tab().append(t.ties).tab().append(h.from).tab().append(h.to).tab().append(clen)
					.tab().append(h.start).tab().append(h.stop).tab().append(h.strand).tab().append(h.trunc).tab().append(h.included ? 1 : 0).tab().append(h.full ? 1 : 0).tab().append(2L*(h.stop-h.start+1)>=clen ? 1 : 0);}
				primary.print(b.nl());add("all", "ALL", t);add(t.clade, "ALL", t);add("all", t.role, t);add(t.clade, t.role, t);
			}
			summary.println("model\tclen\tflank\tclade\trole\ttargets\tbases\treportedHits\tincludedHits\tfullModelHits\ttargetsWithHit\ttargetsWithoutReportedHit\tprimaryFullModel\tprimaryTruncated\ttiedPrimaryTargets");
			for(Group g:groups){require(g.counts[5]+g.counts[6]==g.counts[0], "Every target must be hit-bearing or without a reported hit");b.clear().append(accession).tab().append(clen).tab().append(flank).tab().append(g.clade).tab().append(g.role);for(long v:g.counts){b.tab().append(v);}summary.print(b.nl());}
		}finally{RrnaResourceIO.finish(primary, summary);}
	}
	void add(String clade, String role, Target t){
		assert(t.hits>=t.included && t.included>=t.full):"Full-model included hits are a subset of all reported hits";
		final String key=clade+'\t'+role;Group g=byGroup.get(key);if(g==null){g=new Group(clade, role);byGroup.put(key, g);groups.add(g);}final long[] c=g.counts;
		c[0]++;c[1]+=t.length;c[2]+=t.hits;c[3]+=t.included;c[4]+=t.full;
		if(t.best==null){c[6]++;}else{c[5]++;if(t.best.full){c[7]++;}if(!t.best.trunc.equals("no")){c[8]++;}if(t.ties>1){c[9]++;}}
	}
	/** Order every retained field so equal scores cannot make summaries depend on table row order. */
	static boolean precedes(Hit h, int mf, int mt, int a, int z, char strand, String trunc, boolean included, float score, double evalue){
		assert(h!=null && score==h.score):"The retained-field tie-break only compares hits tied at the highest score";
		if(a!=h.start){return a<h.start;}if(z!=h.stop){return z<h.stop;}if(strand!=h.strand){return strand<h.strand;}
		if(mf!=h.from){return mf<h.from;}if(mt!=h.to){return mt<h.to;}
		final int tr=trunc.compareTo(h.trunc);if(tr!=0){return tr<0;}if(included!=h.included){return !included;}
		final int ev=Double.compare(evalue, h.evalue);if(ev!=0){return ev<0;}
		return Float.compare(score, h.score)<0;// Distinguish signed zero in the printed score too.
	}
	static boolean starts(byte[] row, String prefix){if(row.length<prefix.length()){return false;}for(int i=0; i<prefix.length(); i++){if(row[i]!=prefix.charAt(i)){return false;}}return true;}
	public static String text(byte[] row){return new String(row, StandardCharsets.US_ASCII);}
	static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	static final class Target {
		Target(String i, String c, String a, String r, int n){id=i;clade=c;assembly=a;role=r;length=n;}
		final String id, clade, assembly, role;final int length;long hits, included, full, ties;Hit best;
	}
	static final class Hit {
		Hit(int f, int t, int a, int z, char s, String tr, boolean inc, boolean full_, float score_, double e){from=f;to=t;start=a;stop=z;strand=s;trunc=tr;included=inc;full=full_;score=score_;evalue=e;}
		final int from, to, start, stop;final char strand;final String trunc;final boolean included, full;final float score;final double evalue;
	}
	static final class Group {Group(String c, String r){clade=c;role=r;}final String clade, role;final long[] counts=new long[10];}
	final String accession;final int flank;int clen;long rawHits;
	final ObjectMap<String,Target> byId=new ObjectMap<String,Target>(String.class, Target.class);final ArrayList<Target> ordered=new ArrayList<Target>();
	final ObjectMap<String,Group> byGroup=new ObjectMap<String,Group>(String.class, Group.class);final ArrayList<Group> groups=new ArrayList<Group>();
}
