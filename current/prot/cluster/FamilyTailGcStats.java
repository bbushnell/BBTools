package prot.cluster;

import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.util.*;
import dna.GeneticCode;
import ddl.DDLFamilySearch.StructuralInputs;
import fileIO.*;
import parse.*;
import prot.FamilyGcInputs;
import prok.StopTailCandidateScan;
import shared.Shared;
import stream.*;
import structures.*;

/** Equal-member mean/sample-SD of the exact retained sense-codon tail.
 * Joins frozen input order, memberships, literal extraction metadata and original DNA.
 * One genome is resident; no bulk per-query output or synthesized codons.
 * @author Qiqi
 */
public final class FamilyTailGcStats{
	public static void main(String[] args)throws Exception{
		final PreParser pp=new PreParser(args,FamilyTailGcStats.class,false);
		try{new FamilyTailGcStats(pp.args).run();}finally{Shared.closeStream(pp.outstream);}
	}
	private FamilyTailGcStats(String[] args){
		final HashMap<String,String> o=new HashMap<String,String>();
		for(String a:args){final String[] p=Parse.splitOnFirst(a,'=');require(p.length==2 && o.put(p[0],p[1])==null,"Expected unique flag=value");}
		inputs=take(o,"inputs");members=take(o,"members");familyMap=take(o,"familymap");genomes=take(o,"genomes");readers=take(o,"readers");
		out=Paths.get(take(o,"out"));expectedRecords=Long.parseLong(take(o,"records"));expectedMembers=Long.parseLong(take(o,"assigned"));
		expectedFamilies=Integer.parseInt(take(o,"families"));minAa=Integer.parseInt(take(o,"minaa"));maxAa=Integer.parseInt(take(o,"maxaa"));
		require(o.isEmpty() && (minAa==40 && maxAa==99 || minAa==100 && maxAa==100),"Unknown arguments or unsupported view");
		require(expectedFamilies>0 && expectedMembers>=expectedFamilies && expectedRecords>=expectedMembers,"Invalid declared population");
		require(out.isAbsolute() && !Files.exists(out),"Fresh absolute output directory required");
	}
	private void run()throws Exception{
		loadGenomes();loadReaders();loadFamilies();Files.createDirectory(out);
		final ByteFile in=lines(inputs),mf=lines(members);
		final LineParser1 p=new LineParser1('\t'),m=new LineParser1('\t');
		try{
			header(in,"ordinal\tpath");header(mf,"input_index\tinput_id\tfamily_id\tstatus");
			int ordinal=0;
			for(byte[] row=in.nextLine();row!=null;row=in.nextLine()){
				p.set(row);require(p.terms()==2 && p.parseInt(0)==ordinal++,"Ordered FASTA manifest requires consecutive ordinals");
				final String path=p.parseString(1);tracker.track(Paths.get(path));
				final Streamer aa=stream(path,true);aa.start();
				try{
					for(ListNum<Read> batch=aa.nextList();batch!=null && batch.size()>0;batch=aa.nextList()){
						for(Read r:batch){
							// Clustering membership and extraction metadata use the FASTA ID token;
							// retained tails also carry taxonomic comments after that token.
							int end=0;while(end<r.id.length() && r.id.charAt(end)>32){end++;}
							require(end>0,"Empty protein FASTA identifier");
							if(end<r.id.length()){r.id=r.id.substring(0,end);}
							final byte[] member=mf.nextLine();require(member!=null,"Membership ended before input");
							m.set(member);
							if(m.terms()!=4 || m.parseLong(0)!=records || !m.termEquals(r.id,1)){
								throw new IllegalArgumentException("Membership/input ordinal-ID mismatch at "+records+": "+r.id);
							}
							records++;require(r.bases!=null && r.length()>=minAa && r.length()<=(maxAa==100 ? 200 : 99),"Input length outside retained view");
							final int family=m.parseInt(2);
							if(m.termEquals("assigned",3)){
								final int rank=byId.get(family);require(rank>=0,"Assigned family absent from current map");
								measure(r,rank);assigned++;
							}else{require(family==-1 && m.termEquals("unassigned",3),"Unexpected membership status");}
						}
					}
				}finally{aa.close();}
				require(!aa.errorState(),"Protein input read failed: "+path);
			}
			require(mf.nextLine()==null,"Extra membership rows after input");
		}finally{final boolean a=in.close(),b=mf.close();require(!a && !b,"Input/membership read failed");closeGenome();}
		require(records==expectedRecords && assigned==expectedMembers,"Global input/member coverage differs");
		for(Family f:families){require(f.n==f.expected && f.residues==f.expectedResidues,"Family coverage differs: "+f.id);}
		verifySources();write();verifySources();
		Files.write(out.resolve("PASS"),("FAMILY_TAIL_GC_STATS_PASS records="+records+" assigned="+assigned+" families="+families.size()+
			" missing=0 interval=sense_codons_excluding_stop sample_sd=true\n").getBytes(StandardCharsets.US_ASCII));
		System.out.println(new String(Files.readAllBytes(out.resolve("PASS")),StandardCharsets.US_ASCII));
	}
	private void loadGenomes()throws Exception{
		final ByteFile in=lines(genomes);final LineParser1 p=new LineParser1('\t');
		try{
			header(in,"genome_id\tsource\tfasta\tgff\ttaxid\tphylum_taxid\ttaxonomy_namespace\tsplit\ttable\ttable_assumed");
			for(byte[] row=in.nextLine();row!=null;row=in.nextLine()){
				p.set(row);require(p.terms()==10,"Genome manifest requires ten columns");
				final Genome g=new Genome(p.parseString(0),p.parseString(2));
				require(genomeMap.put(g.id,g)==null,"Duplicate genome identity");
			}
		}finally{require(!in.close(),"Genome manifest failed");}
	}
	private void loadReaders()throws Exception{
		final ByteFile in=lines(readers);final LineParser1 p=new LineParser1('\t');
		try{
			header(in,"genome_id\tsplit\ttails_fasta\tmeta_tsv\tstats_tsv\taudit_tsv");
			for(byte[] row=in.nextLine();row!=null;row=in.nextLine()){
				p.set(row);require(p.terms()==6,"Reader manifest requires six columns");
				final Genome g=genomeMap.get(p.parseString(0));require(g!=null && g.meta==null,"Reader genome missing or repeated");
				g.meta=p.parseString(3);
			}
		}finally{require(!in.close(),"Reader manifest failed");}
	}
	private void loadFamilies()throws Exception{
		final ByteFile in=lines(familyMap);final LineParser1 p=new LineParser1('\t');long n=0;
		try{
			header(in,"stable_family_id\troster_rank\tfamily_id\tassigned\tresidues");
			for(byte[] row=in.nextLine();row!=null;row=in.nextLine()){
				p.set(row);require(p.terms()==5,"Family map requires five columns");
				final int id=p.parseInt(0),rank=p.parseInt(1);final long count=p.parseLong(3),res=p.parseLong(4);
				require(id>=0 && rank==families.size() && byId.get(id)<0 && count>0,"Invalid family identity/order/count");
				byId.put(id,rank);families.add(new Family(id,rank,p.parseString(2),count,res));n=Math.addExact(n,count);
			}
		}finally{require(!in.close(),"Family map failed");}
		require(families.size()==expectedFamilies && n==expectedMembers,"Declared family population differs");
	}
	private void measure(Read r,int rank)throws Exception{
		final int bar=r.id.indexOf('|');require(bar>0,"Missing genome in tail ID");
		if(current==null || current.id.length()!=bar || !r.id.regionMatches(0,current.id,0,bar)){
			closeGenome();current=genomeMap.get(r.id.substring(0,bar));
			require(current!=null && current.meta!=null && !current.visited,"Missing genome metadata or noncontiguous genome input");
			current.visited=true;meta=lines(current.meta);header(meta,META_HEADER);loadDna(current);
			genomesRead++;if(genomesRead%1000==0){System.err.println("TAIL_GC_PROGRESS genomes="+genomesRead+" records="+records+" assigned="+assigned);}
		}
		boolean found=false;
		for(byte[] row=meta.nextLine();row!=null;row=meta.nextLine()){
			metadata.set(row);require(metadata.terms()==24,"Extraction metadata requires24 columns");
			if(metadata.termEquals(r.id,0)){found=true;break;}
		}
		require(found && metadata.termEquals(current.id,1),"Assigned tail missing from exact extraction metadata: "+r.id);
		if(contigName==null || !metadata.termEquals(contigName,2)){
			contigName=metadata.parseString(2);contig=dna.get(decode(contigName));
		}
		require(contig!=null,"Metadata contig absent from source DNA");
		final boolean minus=metadata.termEquals("-",4);
		require(minus || metadata.termEquals("+",4),"Unknown strand");
		final int endpoint=metadata.parseInt(5),stop=minus ? endpoint-1 : endpoint-3;
		final int table=metadata.parseInt(16),original=metadata.parseInt(19),length=metadata.parseInt(20);
		require(length==r.length() && original>=length && metadata.parseInt(22)==0,"Metadata/peptide length or terminal completeness differs");
		final int start=metadata.parseInt(6),end=metadata.parseInt(7),phase=metadata.parseInt(8);
		require(start>=1 && end<=contig.length && phase>=0 && phase<=2 && end-start+1-phase==3*(original+1) &&
			endpoint==(minus ? start : end),"CDS coordinate/phase/length contract differs");
		final GeneticCode code=GeneticCode.forTable(table);
		int codon=GeneticCode.codon(contig,stop);if(minus && codon>=0){codon=rc(codon);}
		require(codon>=0 && code.isStop(codon),"Recorded endpoint is not a real stop");
		final int tailLength=Math.min(maxAa,r.length()),offset=r.length()-tailLength;
		final int gc=StopTailCandidateScan.senseGcCount(contig,stop,minus ? 1 : 0,tailLength);
		for(int i=0;i<tailLength;i++){
			final int pos=minus ? stop+3+3*(tailLength-1-i) : stop-3*tailLength+3*i;
			codon=GeneticCode.codon(contig,pos);if(minus && codon>=0){codon=rc(codon);}
			require(codon>=0 && code.aminoAcid(codon)==r.bases[offset+i],"Actual DNA does not translate to retained peptide: "+r.id);
		}
		families.get(rank).add(gc/(3.0*tailLength),tailLength);
	}
	private void loadDna(Genome genome)throws Exception{
		// Published DNA manifests include aliases into older retained sequence trees.
		// Read and track one resolved regular file; preserve the alias binding as well.
		genome.resolved=Paths.get(genome.fasta).toRealPath();
		tracker.track(genome.resolved);dna.clear();contigName=null;contig=null;
		final Streamer input=stream(genome.resolved.toString(),false);input.start();
		try{
			for(ListNum<Read> batch=input.nextList();batch!=null && batch.size()>0;batch=input.nextList()){
				for(Read r:batch){
					int end=0;while(end<r.id.length() && !Character.isWhitespace(r.id.charAt(end))){end++;}
					require(end>0 && dna.put(r.id.substring(0,end),r.bases)==null,"Duplicate or empty contig identifier");
				}
			}
		}finally{input.close();}
		require(!input.errorState() && !dna.isEmpty(),"DNA FASTA read failed");
	}
	private void verifySources()throws Exception{
		tracker.verify();
		for(Genome g:genomeMap.values()){
			if(g.resolved!=null){require(g.resolved.equals(Paths.get(g.fasta).toRealPath()),"Source DNA alias changed: "+g.id);}
		}
	}
	private void closeGenome(){
		if(meta!=null){require(!meta.close(),"Metadata read failed");meta=null;}
	}
	private void write(){
		final ByteStreamWriter w=new ByteStreamWriter(out.resolve("family_gc.tsv").toString(),false,false,true);w.start();
		final ByteBuilder b=new ByteBuilder();w.println(FamilyGcInputs.TAIL_HEADER_TEXT);
		try{
			for(Family f:families){
				b.clear().append(f.stable).tab().append(f.rank).tab().append(f.id).tab().append(f.expected).tab().append(f.expectedResidues)
					.tab().append(f.n).tab().append(Double.toString(f.mean)).tab().append(Double.toString(f.n>1 ? Math.sqrt(f.m2/(f.n-1)) : 0))
					.tab().append(minAa).tab().append(maxAa).nl();w.print(b);
			}
		}finally{require(!w.poisonAndWait(),"Family GC output failed");}
	}
	private ByteFile lines(String path)throws Exception{tracker.track(Paths.get(path));return ByteFile.makeByteFile1(path,false);}
	private static Streamer stream(String path,boolean amino){
		// The native zero-thread FASTA reader snapshots this flag during construction.
		final boolean old=Shared.AMINO_IN;Shared.AMINO_IN=amino;
		try{return StreamerFactory.makeStreamer(FileFormat.testInput(path,FileFormat.FASTA,null,true,false),null,0,true,-1,false,true,0);}
		finally{Shared.AMINO_IN=old;}
	}
	private static int rc(int c){return ((3^(c&3))<<4)|((3^((c>>>2)&3))<<2)|(3^(c>>>4));}
	/** Reverse the extraction metadata's percent escaping; literal plus remains plus. */
	private static String decode(String s){
		if(s.indexOf('%')<0){return s;}
		final byte[] raw=s.getBytes(StandardCharsets.UTF_8);final ByteBuilder b=new ByteBuilder();
		for(int i=0;i<raw.length;i++){
			if(raw[i]!='%'){b.append(raw[i]);continue;}
			require(i+2<raw.length,"Truncated contig escape");
			final int hi=Character.digit((char)raw[++i],16),lo=Character.digit((char)raw[++i],16);
			require(hi>=0 && lo>=0,"Invalid contig escape");b.append((byte)((hi<<4)|lo));
		}
		return new String(b.toBytes(),StandardCharsets.UTF_8);
	}
	private static void header(ByteFile f,String h){require(Arrays.equals(f.nextLine(),h.getBytes(StandardCharsets.US_ASCII)),"Unexpected header: "+h);}
	private static String take(HashMap<String,String> o,String key){final String v=o.remove(key);require(v!=null && !v.isEmpty(),"Missing "+key);return v;}
	private static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final class Genome{Genome(String i,String f){id=i;fasta=f;}final String id,fasta;String meta;Path resolved;boolean visited;}
	private static final class Family{
		Family(int s,int r,String i,long e,long er){stable=s;rank=r;id=i;expected=e;expectedResidues=er;}
		void add(double x,int length){n++;residues+=length;final double delta=x-mean;mean+=delta/n;m2+=delta*(x-mean);}
		final int stable,rank;final String id;final long expected,expectedResidues;long n,residues;double mean,m2;
	}
	private final String inputs,members,familyMap,genomes,readers;
	private final Path out;private final long expectedRecords,expectedMembers;private final int expectedFamilies,minAa,maxAa;
	private final StructuralInputs tracker=new StructuralInputs();
	private final HashMap<String,Genome> genomeMap=new HashMap<String,Genome>();
	private final HashMap<String,byte[]> dna=new HashMap<String,byte[]>();
	private final ArrayList<Family> families=new ArrayList<Family>();private final IntHashMap byId=new IntHashMap();
	private final LineParser1 metadata=new LineParser1('\t');
	private Genome current;private ByteFile meta;private String contigName;private byte[] contig;
	private long records,assigned,genomesRead;
	private static final String META_HEADER="tail_id\tgenome_id\tcontig\tcontig_ordinal\tstrand\tstop_coordinate\tcds_start\tcds_end\tphase\tnative_id\tnative_aliases\tsource\tgff_source\ttaxid\tphylum_taxid\ttaxonomy_namespace\ttable\tcode_source\ttable_assumed\toriginal_length\ttail_length\tpartial_5prime\tpartial_3prime\tsplit";
}
