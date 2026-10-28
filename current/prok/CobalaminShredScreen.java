package prok;

import java.util.ArrayList;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectMap;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import structures.ByteBuilder;
import structures.IntList;

/** Excludes random cores overlapping full-span, GA-passing clan CM hits.
 * Consumes explicitly selected padded contexts, preserves all raw hit evidence,
 * and distinguishes scanned zero-hit cores from unadjudicated background.
 * @author Raiden
 */
public final class CobalaminShredScreen {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args, CobalaminShredScreen.class, false);
		String manifest=null, cores=null, models=null, out=null, scanned=null;int expected=-1, expectedScanned=-1;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0, "Use flag=value");final String k=arg.substring(0, eq).toLowerCase(Locale.ROOT), v=arg.substring(eq+1);
			if(k.equals("manifest")){manifest=v;}else if(k.equals("in")){cores=v;}else if(k.equals("models")){models=v;}else if(k.equals("outprefix")){out=v;}
			else if(k.equals("scanned")){scanned=v;}else if(k.equals("expectedscanned")){expectedScanned=Integer.parseInt(v);}
			else if(k.equals("expected")){expected=Integer.parseInt(v);}else{throw new IllegalArgumentException("Unknown option: "+k);}
		}
		require(manifest!=null && cores!=null && models!=null && out!=null && scanned!=null && expected>0 && expectedScanned>=0 && expectedScanned<=expected,
			"Bind all cores, coordinate manifest, explicit scanned subset, three CM tables and both denominators");
		Shared.setThreads(1);Shared.TRIM_READ_DESCRIPTION=false;
		try{final CobalaminShredScreen x=new CobalaminShredScreen();x.load(manifest, expected);x.scanned(scanned, expectedScanned);x.scan(models);x.write(cores, out);}finally{Shared.closeStream(pp.outstream);}
	}
	private void load(String path, int expected){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{RrnaResourceIO.header(in, CobalaminShredSampler.HEADER, path);for(byte[] row; (row=in.nextLine())!=null; ){
			p.set(row);require(p.terms()==13, "Retain the full core/context coordinate schema");
			final Member m=new Member(p.parseString(0), row, p.parseInt(10), p.parseInt(11), p.parseInt(12), p.parseInt(9)-p.parseInt(8)+1);
			require(m.start>0 && m.stop>=m.start && m.stop<=m.contextLength && m.stop-m.start+1==m.length && p.parseInt(7)-p.parseInt(6)+1==m.length
				&& p.parseInt(6)-p.parseInt(8)+1==m.start && p.parseInt(7)-p.parseInt(8)+1==m.stop && p.parseInt(8)>0,
				"Every core must map exactly into its padded source context: "+m.id);
			require(byId.put(m.id, m)==null, "Repeated core identity: "+m.id);members.add(m);
		}}finally{require(!in.close(), "Sample manifest read failed");}
		require(members.size()==expected, "All pre-screen candidates must enter the denominator: "+members.size()+" expected="+expected);
	}
	private void scanned(String path, int expected){
		final ByteFile in=RrnaResourceIO.open(path);
		try{for(byte[] row; (row=in.nextLine())!=null; ){
			final String id=new String(row, java.nio.charset.StandardCharsets.US_ASCII);final Member m=byId.get(id);
			require(m!=null && !m.scanned, "Every CM-screened context must bind one unique sampled core: "+id);m.scanned=true;scannedCount++;
		}}finally{require(!in.close(), "Scanned-subset ledger read failed");}
		require(scannedCount==expected, "The declared CM-screened denominator changed");
	}
	private void scan(String path){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');int models=0;
		try{RrnaResourceIO.header(in, MODEL_HEADER, path);for(byte[] row; (row=in.nextLine())!=null; ){
			p.set(row);require(p.terms()==5 && models<3 && p.termEquals(FAMILIES[models], 0), "Screen must include all three clan models in canonical order");
			final String model=p.parseString(0), table=p.parseString(3), hash=p.parseString(4);final int clen=p.parseInt(1);final double ga=Double.parseDouble(p.parseString(2));
			require(clen>0 && Double.isFinite(ga), "Bind explicit CM length and GA");
			if(scannedCount==0){require(table.equals(".") && hash.equals("NA"), "An empty selected subset must explicitly record no CM scan, not a fabricated empty hit table");models++;continue;}
			require(table.startsWith("/") && hash.matches("[0-9a-f]{20}"), "Every nonempty selected subset requires an absolute hashed CM table");
			final RrnaCmTruthBuilder.Criteria criteria=new RrnaCmTruthBuilder.Criteria(model, model, clen, clen, ga);
			final ArrayList<RrnaCmTruthBuilder.Hit> hits=RrnaCmTruthBuilder.readHits(new RrnaCmTruthBuilder.Genome("background", "ALL", table, hash), criteria);
			for(RrnaCmTruthBuilder.Hit hit:hits){
				final Member m=byId.get(hit.seqid);require(m!=null && m.scanned, "CM evidence outside the explicitly scanned context set: "+hit.seqid);
				final int lo=Math.min(hit.seqFrom, hit.seqTo), hi=Math.max(hit.seqFrom, hit.seqTo);
				require(hi<=m.contextLength, "CM coordinates must fit the actual padded context: "+m.id);
				final boolean overlap=lo<=m.stop && hi>=m.start, exclude=overlap && hit.status.equals("QUALIFIED");
				m.hits++;if(overlap){m.overlaps++;}if(exclude){m.excludedBy|=1<<models;}
				evidence.add(new Evidence(m, hit, model, overlap, exclude));
			}
			models++;
		}}finally{require(!in.close(), "CM model/table manifest read failed");}
		require(models==3, "An unscanned model cannot be silently labeled negative");
	}
	private void write(String fasta, String prefix){
		RrnaResourceIO.fresh(prefix, ".cores.fa.gz", ".shreds.tsv", ".disposition.tsv", ".hits.tsv", ".lengths.tsv", ".summary.tsv");
		final ByteStreamWriter kept=RrnaResourceIO.writer(prefix+".shreds.tsv"), disposition=RrnaResourceIO.writer(prefix+".disposition.tsv"), hits=RrnaResourceIO.writer(prefix+".hits.tsv"), lengths=RrnaResourceIO.writer(prefix+".lengths.tsv"), summary=RrnaResourceIO.writer(prefix+".summary.tsv");
		final Euk18sShredSet.Output output=new Euk18sShredSet.Output(prefix+".cores.fa.gz", 0);final ByteBuilder b=new ByteBuilder();final long[] counted=new long[4];
		final IntList before=new IntList(), after=new IntList();
		try{
			RrnaResourceIO.stream(fasta, r->{
				final String id=RrnaResourceIO.first(r.id);final Member m=byId.get(id);
				require(m!=null && !m.seen && r.bases!=null && r.length()==m.length, "Every input core must occur once at its declared length: "+id);m.seen=true;
				counted[0]++;counted[1]+=m.length;before.increment(m.length, 1);
				if(m.excludedBy==0){output.add(new Read(r.bases, null, id, counted[0]-1));counted[2]++;counted[3]+=m.length;after.increment(m.length, 1);}
			});
			require(counted[0]==members.size(), "No missing core can disappear from the screened denominator");
			kept.println(CobalaminShredSampler.HEADER);disposition.println("shred\tlength\tstatus\texcluding_models\traw_hits\tcore_overlap_hits");
			for(Member m:members){
				if(m.excludedBy==0){kept.println(m.row);}
				b.clear().append(m.id).tab().append(m.length).tab().append(m.excludedBy!=0?"EXCLUDED_FULL_CLAN_CM_CORE_OVERLAP":m.scanned?"RETAINED_CM_SCREENED_BACKGROUND":"SOURCE_MASKED_BACKGROUND_NOT_CM_ADJUDICATED").tab();
				int labels=0;if(m.excludedBy==0){b.append('.');}else{for(int i=0; i<3; i++){if((m.excludedBy&(1<<i))!=0){if(labels++>0){b.append(',');}b.append(FAMILIES[i]);}}}
				disposition.print(b.tab().append(m.hits).tab().append(m.overlaps).nl());
			}
			hits.println("shred\tmodel\ttable_line\tmodel_from\tmodel_to\tcontext_from\tcontext_to\tstrand\ttrunc\tbits\tevalue\tcm_status\tcore_overlap\texcludes_core");
			for(Evidence e:evidence){final RrnaCmTruthBuilder.Hit h=e.hit;
				hits.print(b.clear().append(e.member.id).tab().append(e.model).tab().append(h.line).tab().append(h.mdlFrom).tab().append(h.mdlTo).tab().append(h.seqFrom).tab().append(h.seqTo)
					.tab().append(h.strand).tab().append(h.trunc).tab().append(h.score).tab().append(h.evalue).tab().append(h.status).tab().append(e.overlap).tab().append(e.exclude).nl());
			}
			lengths.println("length\tcandidates\tretained\texcluded");for(int i=0; i<before.size; i++){
				if(before.get(i)>0){lengths.print(b.clear().append(i).tab().append(before.get(i)).tab().append(after.get(i)).tab().append(before.get(i)-after.get(i)).nl());}
			}
			summary.println("candidates\tcandidate_bases\tretained\tretained_bases\texcluded\traw_cm_hits\tscanned\tunadjudicated\tinterpretation");
			summary.print(b.clear().append(counted[0]).tab().append(counted[1]).tab().append(counted[2]).tab().append(counted[3]).tab().append(counted[0]-counted[2]).tab().append(evidence.size())
				.tab().append(scannedCount).tab().append(members.size()-scannedCount).tab().append(scannedCount==0?"CM_NOT_RUN_EMPTY_SELECTION":"TARGETED_THREE_CM_FULL_SPAN_GA_NO_TRUNC_CORE_OVERLAP_SCREEN")
				.append("_NOT_PROVEN_BIOLOGICAL_NEGATIVES\n"));
		}finally{final boolean failed=output.close();RrnaResourceIO.finish(kept, disposition, hits, lengths, summary);require(!failed, "Screened FASTA output failed");}
		System.err.println("COBALAMIN_SHRED_SCREEN_PASS candidates="+counted[0]+" retained="+counted[2]+" excluded="+(counted[0]-counted[2]));
	}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final class Member {
		Member(String i, byte[] r, int a, int z, int n, int c){id=i;row=r;start=a;stop=z;length=n;contextLength=c;}
		final String id;final byte[] row;final int start, stop, length, contextLength;int hits, overlaps, excludedBy;boolean seen, scanned;
	}
	private static final class Evidence {
		Evidence(Member m, RrnaCmTruthBuilder.Hit h, String f, boolean o, boolean x){member=m;hit=h;model=f;overlap=o;exclude=x;}
		final Member member;final RrnaCmTruthBuilder.Hit hit;final String model;final boolean overlap, exclude;
	}
	private final ArrayList<Member> members=new ArrayList<Member>();private final ObjectMap<String,Member> byId=new ObjectMap<String,Member>(String.class, Member.class);
	private final ArrayList<Evidence> evidence=new ArrayList<Evidence>();
	private int scannedCount;
	private static final String[] FAMILIES={"RF00174", "RF01482", "RF01689"};
	public static final String MODEL_HEADER="model\tclen\tga_bits\ttblout\ttblout_sha80";
}
