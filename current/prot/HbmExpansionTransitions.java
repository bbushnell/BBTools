package prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.List;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.LongHashMap;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Streams matched shipping/restricted/expanded calls into explicit transition counts. @author Keqing */
public final class HbmExpansionTransitions {

	public static void main(String[] args){
		try{
			Shared.setThreads(1);
			final PreParser pp=new PreParser(args, HbmExpansionTransitions.class, false);
			try{run(HmmComparisonData.options(pp.args));}finally{Shared.closeStream(pp.outstream);}
		}catch(Throwable failure){failure.printStackTrace(); System.exit(1);}
	}

	private static void run(HashMap<String,String> o) throws Exception{
		final String[] keys={"identity", "shipping", "restricted", "expanded", "map"};
		for(String key:o.keySet()){
			boolean allowed=key.equals("out") || key.equals("expected") || key.equals("oldfamilies");
			for(String k:keys){allowed|=key.equals(k) || key.equals(k+"sha80");}
			require(allowed, "Unknown transition argument: "+key);
		}
		for(String key:keys){HbmProfilePilot.requireHash(HmmComparisonData.required(o,key), HmmComparisonData.required(o,key+"sha80"));}
		final int expected=parse.Parse.parseIntKMG(HmmComparisonData.required(o,"expected"));
		final int old=parse.Parse.parseIntKMG(HmmComparisonData.required(o,"oldfamilies"));
		final Families families=new Families(o.get("identity"));
		require(expected>0 && old>0 && old<families.names.length, "Positive query count and a strict old-family prefix are required");
		final Path out=Paths.get(HmmComparisonData.required(o,"out"));
		require(!Files.exists(out), "Transition output must be fresh: "+out);
		final Cursor shipping=new Cursor(o.get("shipping"), SHIPPING_HEADER);
		final Cursor restricted=new Cursor(o.get("restricted"), CALL_HEADER);
		final Cursor expanded=new Cursor(o.get("expanded"), CALL_HEADER);
		final Cursor mapping=new Cursor(o.get("map"), "query\tsample\toriginal_id\tlength");
		Files.createDirectory(out);
		final ByteStreamWriter changes=HmmComparisonData.writer(out.resolve("changes.tsv").toString());
		changes.println("query\tsample\toriginal_id\tshipping_family\trestricted_family\texpanded_family\tshipping_to_full\tcontrol_to_full\tshipping_reason\tshipping_identity\tshipping_overlap\trestricted_score64\trestricted_identity\trestricted_mutual_coverage\tfull_score64\tfull_identity\tfull_mutual_coverage");
		final String[] groups=groups(); final long[][][] counts=new long[3][groups.length][CATEGORIES.length];
		final ByteBuilder line=new ByteBuilder(); int n=0;
		try{
			while(mapping.next()){
				require(shipping.next() && restricted.next() && expanded.next(), "Assignment stream ended before the query map at row "+n);
				query(mapping.row,n); query(shipping.row,n); query(restricted.row,n); query(expanded.row,n);
				require(mapping.row.terms()==4 && shipping.row.terms()==10, "Unexpected map/shipping field count at query "+n);
				final int group=group(mapping.row,groups), domain=group<53 ? 1 : 2;
				final int a=shipping.row.termEquals("NA",1) ? -1 : families.find(shipping.row,1);
				require(a<old, "Shipping assignment is outside the original roster at query "+n);
				final int b=assigned(restricted.row,families,old), c=assigned(expanded.row,families,families.names.length);
				final int length=mapping.row.parseInt(3);
				require(length>0 && restricted.row.parseInt(10)==length && expanded.row.parseInt(10)==length, "Matched query lengths differ at query "+n);
				final int ab=category(a,b,old), ac=category(a,c,old), bc=category(b,c,old);
				add(counts[0],ac,group,domain); add(counts[1],bc,group,domain); add(counts[2],ab,group,domain);
				if(a!=b || b!=c){
					line.clear(); field(line,mapping.row,0); line.tab().append(groups[group]).tab(); field(line,mapping.row,2);
					line.tab().append(a<0 ? "NA" : families.names[a]).tab().append(b<0 ? "NA" : families.names[b])
						.tab().append(c<0 ? "NA" : families.names[c]).tab().append(CATEGORIES[ac]).tab().append(CATEGORIES[bc]).tab();
					field(line,shipping.row,8); line.tab(); field(line,shipping.row,3); line.tab(); field(line,shipping.row,4);
				for(int f:METRIC_FIELDS){line.tab(); field(line,restricted.row,f);}
				for(int f:METRIC_FIELDS){line.tab(); field(line,expanded.row,f);}
					changes.print(line.nl());
				}
				n++;
			}
			require(n==expected && !shipping.next() && !restricted.next() && !expanded.next(), "Assignment/map cardinality differs from expected="+expected+" actual="+n);
		}finally{shipping.close(); restricted.close(); expanded.close(); mapping.close(); HmmComparisonData.close(changes);}
		final ByteStreamWriter summary=HmmComparisonData.writer(out.resolve("summary.tsv").toString());
		line.clear().append("comparison\tgroup\tqueries"); for(String c:CATEGORIES){line.tab().append(c);} summary.print(line.nl());
		for(int p=0;p<3;p++){
			for(int g=0;g<groups.length;g++){
				long total=0; for(long value:counts[p][g]){total+=value;}
				if(g==0){require(total==expected, "Transition categories do not conserve every query");}
				line.clear().append(PAIRS[p]).tab().append(groups[g]).tab().append(total);
				for(long value:counts[p][g]){line.tab().append(value);} summary.print(line.nl());
			}
		}
		HmmComparisonData.close(summary);
		for(String key:keys){HbmProfilePilot.requireHash(o.get(key),o.get(key+"sha80"));}
		final ByteStreamWriter pass=HmmComparisonData.writer(out.resolve("PASS").toString());
		pass.println("HBM_EXPANSION_TRANSITIONS_PASS queries="+n+" comparisons=3"); HmmComparisonData.close(pass);
	}

	/** Disjoint outcomes; a new-family switch requires an earlier accepted assignment. */
	static int category(int before,int after,int old){
		require(before>=-1 && after>=-1 && old>0,"Invalid transition identity");
		if(before<0){return after<0 ? 0 : after<old ? 2 : 3;}
		return after==before ? 1 : after<0 ? 4 : after<old ? 5 : 6;
	}
	private static void add(long[][] counts,int category,int group,int domain){
		assert(group>=3 && domain>=1 && domain<=2) : "Groups distinguish all, domain and individual genomes";
		counts[0][category]++; counts[domain][category]++; counts[group][category]++;
	}
	private static String[] groups(){
		final String[] names=new String[103]; names[0]="all"; names[1]="bacteria"; names[2]="archaea";
		for(int i=0;i<50;i++){names[i+3]="bacteria_"+i; names[i+53]="archaea_"+i;} return names;
	}
	private static int group(LineParser1 row,String[] names){
		final boolean bacteria=row.termStartsWith("bacteria_",1);
		require(bacteria || row.termStartsWith("archaea_",1),"Unknown source genome group");
		final int number=row.parseInt(1,bacteria ? 9 : 8), group=number+(bacteria ? 3 : 53);
		require(number>=0 && number<50 && row.termEquals(names[group],1),"Source genome is outside the recorded100-genome pool");
		return group;
	}
	private static void query(LineParser1 row,int expected){
		int digits=1; for(int value=expected;value>=10;value/=10){digits++;}
		require(row.termStartsWith("comparison_g",0) && row.length(0)==12+digits && row.parseInt(0,12)==expected,
			"Missing, reordered or noncanonical query identity at ordinal "+expected);
	}
	private static int assigned(LineParser1 row,Families families,int limit){
		require(row.terms()==13,"Experimental assignment must have13 fields");
		final int rank=row.parseInt(2);
		if(row.termEquals("NO_ELIGIBLE_FAMILY",1)){
			require(rank==-1 && row.parseInt(3)==-1 && row.termEquals("NA",4),"Rejected query carries a family identity"); return -1;
		}
		require(row.termEquals("ASSIGNED",1) && rank>=0 && rank<limit,"Unexpected assignment status or family outside candidate prefix");
		require(row.parseInt(3)==families.permanent[rank] && row.termEquals(families.names[rank],4),"Dense/permanent/representative family IDs differ");
		return rank;
	}
	private static void field(ByteBuilder out,LineParser1 row,int field){
		final int len=row.setBounds(field); out.append(row.line(),row.a(),len);
	}
	/** A lookup accelerator only; every hit is compared against the full retained identifier. */
	private static long hash(byte[] bytes,int from,int to){
		assert(from>=0 && to<=bytes.length && from<to) : "Family hash needs a nonempty bounded field";
		long hash=0xcbf29ce484222325L; for(int i=from;i<to;i++){hash=(hash^(bytes[i]&255))*0x100000001b3L;} return hash;
	}
	private static final class Families {
		Families(String path){
			final List<String[]> rows=HmmComparisonData.rows(path); require(rows.size()>1,"Empty family identity table");
			names=new String[rows.size()-1]; permanent=new int[names.length]; ranks=new LongHashMap(names.length*2);
			require(String.join("\t",rows.get(0)).equals("active_index\tfamily_id\tsource_active_index\trep_id\tmembers\traw_residues"),"Unsupported family identity schema");
			for(int i=0;i<names.length;i++){
				final String[] row=rows.get(i+1); require(row.length==6 && Integer.parseInt(row[0])==i,"Non-dense family identity table");
				names[i]=row[3]; permanent[i]=Integer.parseInt(row[1]); final byte[] bytes=names[i].getBytes(StandardCharsets.US_ASCII);
				final long key=hash(bytes,0,bytes.length); require(!ranks.contains(key),"Duplicate family ID or lookup-hash collision"); ranks.put(key,i);
			}
		}
		int find(LineParser1 row,int field){
			row.setBounds(field); final int rank=ranks.get(hash(row.line(),row.a(),row.b()));
			require(rank>=0 && row.termEquals(names[rank],field),"Shipping family is absent from the final identity table"); return rank;
		}
		final String[] names; final int[] permanent; final LongHashMap ranks;
	}
	private static final class Cursor {
		Cursor(String path,String header){
			file=ByteFile.makeByteFile(path,false); final byte[] first=file.nextLine();
			require(first!=null && new String(first,StandardCharsets.UTF_8).equals(header),"Unexpected assignment/map header: "+path);
		}
		boolean next(){final byte[] line=file.nextLine(); if(line==null){return false;} row.set(line); return true;}
		void close(){require(!file.close(),"Transition input I/O failed");}
		final ByteFile file; final LineParser1 row=new LineParser1('\t');
	}
	private static void require(boolean ok,String message){if(!ok){throw new IllegalArgumentException(message);}}
	private static final String[] CATEGORIES={"both_rejected","unchanged","gained_old","gained_new","lost","switched_old","switched_new"};
	private static final String[] PAIRS={"shipping_to_full","restricted_to_full","shipping_to_restricted"};
	private static final int[] METRIC_FIELDS={5,6,9};
	private static final String SHIPPING_HEADER="query\tfamily\tscore\tidentity\taligned_overlap\traw\thbm_raw\thbm_relative\treason\tcandidates";
	private static final String CALL_HEADER="query\tstatus\tactive_index\tfamily_id\trep_id\tscore64\tidentity_pct\tcoverage_q\tcoverage_core\tmutual_coverage\tquery_length\tref_start\tref_end";
}
