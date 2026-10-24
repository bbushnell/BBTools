package prok;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import map.ObjectIntMap;
import map.ObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import stream.Writer;
import stream.WriterFactory;
import structures.ByteBuilder;
import structures.DoubleList;

/** Descriptive saved-evidence tables only: never changes membership or runs searches.
 * All distributions count annotation copies, not unique sequences. Scores absent
 * from cut_ga scans remain NA, never zero or an inferred subthreshold score.
 * @author Ganyu */
public final class NcrnaMitoCharacterization {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args,NcrnaMitoCharacterization.class,false);
		final Parser parser=new Parser();final NcrnaMitoCharacterization x=new NcrnaMitoCharacterization();String verified=null;
		try{
			for(String arg:pp.args){final String[] kv=Parse.splitOnFirst(arg,'=');final String k=kv[0].toLowerCase(Locale.ROOT),v=kv.length>1?kv[1]:null;
				if(k.equals("outprefix")){x.prefix=v;}else if(k.equals("verified")){verified=v;}else if(k.equals("source")){x.source=v;}
				else if(k.equals("expected")){x.expected=Parse.parseIntKMG(v);}else if(k.equals("expectedcm")){x.expectedCm=Parse.parseIntKMG(v);}
				else if(k.equals("expectednohit")){x.expectedNoHit=Parse.parseIntKMG(v);}else if(k.equals("expected16")){x.expected16=Parse.parseIntKMG(v);}
				else if(k.equals("expected18")){x.expected18=Parse.parseIntKMG(v);}else if(k.equals("lengthbins")){x.bins=parseBins(v);}
				else if(k.equals("expectedssu")){x.expectedSsu=Parse.parseIntKMG(v);}
				else if(!parser.parse(arg,k,v)){throw new IllegalArgumentException("Unknown flag: "+arg);}
			}
			require(parser.in1!=null && x.source!=null && x.prefix!=null && verified!=null && !parser.append && parser.maxReads<0,"Require full verified audit/source FASTA, fresh outprefix and no truncated/append input");
			require(x.expected>0 && x.expectedCm>=0 && x.expectedNoHit>0 && x.expected16>0 && x.expected18>0 && x.expectedSsu>=x.expected16+x.expected18,"Declared population and cohort counts, including the unsplit SSU union, are required");
			final ObjectMap<String,String> pins=EukPanelGffValidation.hashes(verified);
			require(pins.get(parser.in1)!=null && pins.get(x.source)!=null,"The producer must bind both accepted audit and original source FASTA");
			Shared.TRIM_READ_DESCRIPTION=false;ReadWrite.setZipThreads(Math.min(4,Shared.threads()));
			RrnaResourceIO.fresh(x.prefix,".nohit.tsv",".margins.tsv",".nohit.fa.gz",".ssu.fa.gz");x.read(parser.in1);x.loadSequences();x.write();
			x.writeFasta(x.members,".nohit.fa.gz",x.expectedNoHit);x.writeFasta(x.ssuMembers,".ssu.fa.gz",x.expectedSsu);
			pp.outstream.println("MITO_CHARACTERIZATION_PASS annotations="+x.annotations+" cm_only="+x.cmOnly+" nohit="+x.noHit+" m16S="+x.margins[0][0].n()+" m18S="+x.margins[1][0].n()+" exported_nohit="+x.members.size()+" exported_ssu="+x.ssuMembers.size()+" weighting=annotation_copy below_ga=NA no_search=true locked_cohorts_unchanged=true");
		}finally{Shared.closeStream(pp.outstream);}
	}
	private void read(String path){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{
			RrnaResourceIO.header(in,HEADER,path);
			for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==30,"Construction audit schema changed");
				if(p.termEquals("CM_ONLY_EXACT_GEOMETRY",0)){require(p.termEquals("AUDIT_ONLY",4),"CM-only rows cannot enter annotation cohorts");cmOnly++;continue;}
				require(p.termEquals("ANNOTATION",0),"Unknown audit row kind");final String id=p.parseString(1);
				require(!id.isEmpty() && sourceLengths.get(id)<0,"Every source annotation must appear once: "+id);annotations++;
				final int length=p.parseInt(6);require(length>0,"Annotation lengths must remain positive");
				sourceLengths.put(id,length);
				final double a=p.parseDouble(16),b=p.parseDouble(17);require(Double.isFinite(a) && Double.isFinite(b) && a>=0 && a<=1 && b>=0 && b<=1,"Stored universal identities must be finite fractions");
				final int cohort=p.termEquals("m16S",4)?0:p.termEquals("m18S",4)?1:-1;
				if(cohort>=0){
					require(p.termEquals("SSU_SUPPORT_REVIEW",5) && (cohort==0?a>b:a<b),"Locked construction membership must match SSU-only support and signed universal preference");
					int stratum=-1;for(int s=1;s<STRATA.length;s++){if(p.termEquals(STRATA[s],23)){stratum=s;break;}}
					require(stratum>=1,"Unknown CM-concordance stratum in a selected cohort");
					margins[cohort][0].add(a,b,length);margins[cohort][stratum].add(a,b,length);
				}else{require(p.termEquals("UNRESOLVED_TIE",4) || p.termEquals("UNRESOLVED_CONFLICT",4) || p.termEquals("UNRESOLVED_UNSUPPORTED",4) || p.termEquals("LSU_ONLY",4),"Unknown construction class");}
				if(cohort>=0 || p.termEquals("UNRESOLVED_TIE",4)){
					require(p.termEquals("SSU_SUPPORT_REVIEW",5) && (cohort>=0 || a==b),"Unsplit union contains every SSU-only locus including exact preference ties, never conflicts");
					ssuMembers.add(member(p,id));
				}
				if(p.termEquals("NO_GA_HIT_EITHER_STRAND",29)){
					require(p.termEquals("UNRESOLVED_UNSUPPORTED",4) && p.termEquals("NO_SAME_STRAND_GA_HIT_UNRESOLVED",5) && p.termEquals("NA",18) && p.termEquals("NA",19),"No-hit characterization must not acquire a best reported SSU hit or change the selected source subset");
					final String product=p.parseString(9);final int bin=bin(length,bins);final String key=product+'\t'+bin;
					Group g=byGroup.get(key);if(g==null){g=new Group(product,bin);byGroup.put(key,g);groups.add(g);}g.stats.add(a,b,length);noHit++;
					members.add(member(p,id));
				}
			}
		}finally{require(!in.close(),"Audit input read failed");}
		require(annotations==expected && cmOnly==expectedCm && noHit==expectedNoHit && margins[0][0].n()==expected16 && margins[1][0].n()==expected18 && ssuMembers.size()==expectedSsu,"Characterization must conserve all declared source/cohort/no-hit/unsplit denominators");
		for(int c=0;c<2;c++){int n=0;for(int s=1;s<STRATA.length;s++){n+=margins[c][s].n();}require(n==margins[c][0].n(),"Concordance strata must partition each complete cohort");}
	}
	private Member member(LineParser1 p,String id){
		final Member m=new Member();m.id=id;m.unique=p.parseString(12);m.taxid=p.parseInt(3);m.partial=flag(p,7);m.pseudo=flag(p,8);
		require(m.taxid>0 && !m.unique.isEmpty() && byMember.get(id)==null,"Exported cohorts need positive host taxids, original unique IDs and disjoint annotation membership");byMember.put(id,m);return m;
	}
	private void loadSequences(){
		assert(members.size()==expectedNoHit):"Only the complete no-GA audit subset is exported as the distinct third construction cohort";
		final int[] count={0};RrnaResourceIO.stream(source,r->{final String id=RrnaResourceIO.first(r.id);final int length=sourceLengths.get(id);
			require(length>0 && r.bases!=null && r.length()==length,"Source FASTA must contain every original annotation exactly once with unchanged length: "+id);
			sourceLengths.put(id,-length);count[0]++;final Member m=byMember.get(id);if(m!=null){m.bases=r.bases.clone();}
		});require(count[0]==expected,"Complete source FASTA count must equal the accepted annotation census");
		for(Member m:members){require(m.bases!=null,"Missing no-GA FASTA member: "+m.id);}
		for(Member m:ssuMembers){require(m.bases!=null,"Missing unsplit SSU FASTA member: "+m.id);}
	}
	private void writeFasta(ArrayList<Member> selected,String suffix,int expectedCount){
		assert(selected.size()==expectedCount):"Export cannot change the declared source cohort count";
		final Writer writer=WriterFactory.makeWriter(FileFormat.testOutput(prefix+suffix,FileFormat.FASTA,null,true,false,false,true),null,1);writer.start();
		ArrayList<Read> batch=new ArrayList<Read>(256);long batchId=0,count=0;
		try{for(Member m:selected){final String id=m.id+" unique="+m.unique+" tid|"+m.taxid+" partial="+m.partial+" pseudo="+m.pseudo;
			batch.add(new Read(m.bases,null,id,count++,false));if(batch.size()==256){writer.add(batch,batchId++);batch=new ArrayList<Read>(256);}
		}if(!batch.isEmpty()){writer.add(batch,batchId);}}finally{require(!writer.poisonAndWait(),"Third-cohort FASTA output failed");}
		require(writer.readsWritten()==expectedCount,"Every selected source copy must be retained without dropping ties or flagged members");
	}
	private void write(){
		assert(noHit==expectedNoHit):"Only the complete declared no-hit subset can produce a characterization table";
		groups.sort(Comparator.comparing((Group g)->g.product).thenComparingInt(g->g.bin));
		final ByteStreamWriter no=RrnaResourceIO.writer(prefix+".nohit.tsv"),margin=RrnaResourceIO.writer(prefix+".margins.tsv");final ByteBuilder b=new ByteBuilder();
		try{
			no.println("source_product\tlength_bin_inclusive\tannotations\tnohit_denominator\tmin_length\tmax_length\tbest_below_ga_ssu_bits\tbest_below_ga_lsu_bits\tscore_missing_reason\tavailable_cm_evidence\tmean_identity_16S\tp05_identity_16S\tp50_identity_16S\tp95_identity_16S\tmean_identity_18S\tp05_identity_18S\tp50_identity_18S\tp95_identity_18S\tmean_16S_minus_18S\tweighting\tquantile_method");
			int total=0;for(Group g:groups){final Stats s=g.stats;s.sort();total+=s.n();
				b.clear().append(g.product).tab().append(binName(g.bin,bins)).tab().append(s.n()).tab().append(noHit).tab().append(s.minLength).tab().append(s.maxLength)
					.append("\tNA\tNA\tNOT_RETAINED_BY_CUT_GA\tNO_REPORTED_GA_HIT_EITHER_STRAND");
				appendIdentity(b,s.first,s.sumFirst);appendIdentity(b,s.second,s.sumSecond);
				no.print(b.tab().append(s.sumMargin/s.n(),9).append("\tANNOTATION_COPY\tNEAREST_RANK_CEIL_P_TIMES_N").nl());
			}require(total==noHit,"Product/length bins must conserve the exact no-hit annotation subset");
			margin.println("construction\tcm_concordance\tannotations\tcohort_denominator\tmean_identity_16S\tmean_identity_18S\tmin_16S_minus_18S\tp05_margin\tp25_margin\tp50_margin\tp75_margin\tp95_margin\tmax_margin\tmean_margin\tabs_margin_le_0.01\tabs_margin_le_0.02\tabs_margin_le_0.05\tweighting\tquantile_method");
			for(int c=0;c<2;c++){for(int t=0;t<STRATA.length;t++){final Stats s=margins[c][t];s.sort();b.clear().append(c==0?"m16S":"m18S").tab().append(STRATA[t]).tab().append(s.n()).tab().append(margins[c][0].n());
				if(s.n()==0){for(int i=0;i<10;i++){b.append("\tNA");}}
				else{
					b.tab().append(s.sumFirst/s.n(),9).tab().append(s.sumSecond/s.n(),9).tab().append(s.margin.array[0],9);
					for(double q:QUANTILES){b.tab().append(quantile(s.margin,q),9);}
					b.tab().append(s.margin.array[s.n()-1],9).tab().append(s.sumMargin/s.n(),9);
				}
				margin.print(b.tab().append(s.close01).tab().append(s.close02).tab().append(s.close05).append("\tANNOTATION_COPY\tNEAREST_RANK_CEIL_P_TIMES_N").nl());
			}}
		}finally{RrnaResourceIO.finish(no,margin);}
	}
	private static void appendIdentity(ByteBuilder b,DoubleList values,double sum){
		assert(values.size>0):"Only occupied product/length bins enter identity summaries";
		b.tab().append(sum/values.size,9).tab().append(quantile(values,.05),9).tab().append(quantile(values,.5),9).tab().append(quantile(values,.95),9);
	}
	public static double quantile(DoubleList sorted,double q){
		require(sorted!=null && sorted.size>0 && q>0 && q<=1,"Nearest-rank quantiles require a nonempty sorted sample and 0<p<=1");
		return sorted.array[Math.max(0,(int)Math.ceil(q*sorted.size)-1)];
	}
	public static int[] parseBins(String value){
		require(value!=null && !value.isEmpty(),"Length bins must be explicit positive ascending boundaries");final String[] tokens=value.split(",");final int[] result=new int[tokens.length];
		for(int i=0;i<result.length;i++){result[i]=Integer.parseInt(tokens[i]);require(result[i]>1 && (i==0 || result[i]>result[i-1]),"Length-bin boundaries must be strictly ascending and greater than1");}return result;
	}
	public static int bin(int length,int[] limits){require(length>0 && limits!=null,"Positive source length and defined descriptive bins required");int b=0;while(b<limits.length && length>=limits[b]){b++;}return b;}
	private static String binName(int b,int[] limits){assert(b>=0 && b<=limits.length):"Bin index comes from the validated ordered boundary list";return (b==0?1:limits[b-1])+"-"+(b==limits.length?"INF":Integer.toString(limits[b]-1));}
	private static void require(boolean value,String message){if(!value){throw new IllegalArgumentException(message);}}
	private static boolean flag(LineParser1 p,int field){require(p.termEquals("true",field) || p.termEquals("false",field),"Source flags must remain literal booleans");return p.termEquals("true",field);}
	private static final class Stats {
		void add(double a,double b,int length){assert(length>0 && Double.isFinite(a) && Double.isFinite(b)):"Only validated audit identities and lengths enter statistics";
			final double d=a-b;first.add(a);second.add(b);margin.add(d);sumFirst+=a;sumSecond+=b;sumMargin+=d;minLength=Math.min(minLength,length);maxLength=Math.max(maxLength,length);
			// Input identities have nine decimal places. This epsilon only absorbs binary subtraction rounding at a declared threshold.
			if(Math.abs(d)<=.01+1e-12){close01++;}if(Math.abs(d)<=.02+1e-12){close02++;}if(Math.abs(d)<=.05+1e-12){close05++;}
		}
		int n(){return margin.size;}void sort(){assert(first.size==second.size && first.size==margin.size):"The three distributions share exactly the same annotation denominator";first.sort();second.sort();margin.sort();}
		final DoubleList first=new DoubleList(16),second=new DoubleList(16),margin=new DoubleList(16);
		double sumFirst,sumSecond,sumMargin;int minLength=Integer.MAX_VALUE,maxLength,close01,close02,close05;
	}
	private static final class Group {Group(String p,int b){product=p;bin=b;}final String product;final int bin;final Stats stats=new Stats();}
	private static final class Member {String id,unique;int taxid;boolean partial,pseudo;byte[] bases;}
	private NcrnaMitoCharacterization(){for(int c=0;c<2;c++){for(int s=0;s<STRATA.length;s++){margins[c][s]=new Stats();}}}
	private String prefix,source;private int expected=37716,expectedCm=718,expectedNoHit=7050,expected16=4831,expected18=6586,expectedSsu=11420,annotations,cmOnly,noHit;
	private int[] bins={200,500,1000,1500,2000,3000,5000};private final Stats[][] margins=new Stats[2][STRATA.length];
	private final ObjectMap<String,Group> byGroup=new ObjectMap<String,Group>(String.class,Group.class);private final ArrayList<Group> groups=new ArrayList<Group>();
	private final ObjectIntMap<String> sourceLengths=new ObjectIntMap<String>(String.class);
	private final ObjectMap<String,Member> byMember=new ObjectMap<String,Member>(String.class,Member.class);private final ArrayList<Member> members=new ArrayList<Member>();
	private final ArrayList<Member> ssuMembers=new ArrayList<Member>();
	private static final double[] QUANTILES={.05,.25,.5,.75,.95};
	private static final String[] STRATA={"ALL","CONCORDANT_NOMINAL_MODEL_CLASS","DISCORDANT_NOMINAL_MODEL_CLASS","SPECIALIZED_MITO_MODEL","TIED_BEST_MODELS"};
	public static final String HEADER="kind\tid\taccession\ttaxid\tconstruction\tcm_status\tlength\tpartial\tpseudo\tproduct\tsource_segments\tsource_crosses_origin\tunique_sequence_id\tsource_copies\tconstruction_copies\tuniversal_preference\tquantum_16S\tquantum_18S\tbest_ssu_cm\tbest_ssu_bits\tmodel_coverage\ttruncation\tbest_tied_models\tcm_concordance\tsame_strand_overlap_bp\tsource_overlap_fraction\thit_overlap_fraction\thit_segments\texact_genomic_evidence_members\tunsupported_description";
}
