package prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Paired export of pinned raw calibration tables; never recalculates biological scores.
 * Raw tables are immutable. Retained integers become guaranteed minima after lower rounding.
 * @author Sayu
 */
public final class FamilyGateExport{

	private FamilyGateExport(){}

	public static void main(String[] args){
		try{
			final PreParser pp=new PreParser(args,FamilyGateExport.class,false);
			try{run(pp.args);}finally{Shared.closeStream(pp.outstream);}
		}catch(Throwable failure){failure.printStackTrace(); System.exit(1);}
	}

	public static void run(final String[] args) throws Exception{
		final HashMap<String,String> options=new HashMap<String,String>();
		for(String arg : args){
			final String[] pair=Parse.splitOnFirst(arg,'=');
			require(pair.length==2,"Require flag=value arguments");
			final String key=pair[0].toLowerCase(Locale.ROOT);
			require(Arrays.asList("lower","lowersha80","upper","uppersha80","out","maxfamilies").contains(key) &&
				options.put(key,pair[1])==null,"Unknown or duplicate option: "+key);
		}
		final int maxFamilies=Integer.parseInt(options.containsKey("maxfamilies") ? options.get("maxfamilies") : "1000000");
		require(maxFamilies>0,"maxfamilies must be positive");
		final Path lower=Paths.get(required(options,"lower")), upper=Paths.get(required(options,"upper"));
		final String lowerPin=DigestSuffix.requireSuffix(required(options,"lowersha80"),"raw lower"),
			upperPin=DigestSuffix.requireSuffix(required(options,"uppersha80"),"raw upper");
		pin(lower,lowerPin); pin(upper,upperPin);
		require(!Files.isSameFile(lower,upper),"Lower and upper inputs must be distinct");
		final Path out=Paths.get(required(options,"out"));
		require(out.isAbsolute() && out.normalize().equals(out) && !Files.exists(out) && !Files.isSymbolicLink(out) &&
			out.getParent()!=null && Files.isDirectory(out.getParent()) && out.getParent().toRealPath().equals(out.getParent()),
			"Require a fresh absolute output directory with an existing real parent");
		final ArrayList<Row> rows=new ArrayList<Row>(); final HashMap<String,Row> byId=new HashMap<String,Row>();
		readLower(lower,rows,byId,maxFamilies);
		Collections.sort(rows,new Comparator<Row>(){public int compare(Row a,Row b){return Integer.compare(a.rank,b.rank);}});
		int previous=-1;
		for(Row row : rows){require(row.rank>previous,"Duplicate raw lower rank: "+row.rank); previous=row.rank;}
		readUpper(upper,byId);
		long members=0; int zeros=0; final int[] disabled=new int[3];
		for(Row row : rows){
			members=Math.addExact(members,row.members);
			if(row.members==0){zeros++; continue;}
			require(row.hasUpper,"Missing paired upper family: "+row.id);
			for(int i=0; i<3; i++){
				row.lowerEncoded[i]=FamilyGateCodec.lower(row.lower[i]);
				row.upperEncoded[i]=i==0 ? FamilyGateCodec.cosineAccept(row.lower[i],row.upper[i]) :
					i==2 ? FamilyGateCodec.finalBlosum(row.lower[i],row.upper[i]) :
					FamilyGateCodec.upper(row.lower[i],row.upper[i]);
				if(i==1 && FamilyGateCodec.disabled(row.upperEncoded[i])){disabled[i]++;}
			}
		}
		// Check the actual pinned bytes again after parsing, before publishing a new directory.
		pin(lower,lowerPin); pin(upper,upperPin);
		Files.createDirectory(out);
		writeTables(out,rows);
		write(out.resolve("recipe.txt"),"schema="+FamilyGateCodec.SCHEMA+"\nraw_lower="+lower+"\nraw_upper="+upper+
			"\nraw_lower_sha80="+lowerPin+"\nraw_upper_sha80="+upperPin+
			"\nlower_rounding=toward_negative_infinity\nupper_rounding=toward_positive_infinity\n"+
			"identity_final=round_up_float32(max(raw_lower,1.1*raw_negative_upper,raw_negative_upper+2,30));accept_strictly_greater\n"+
			"cosine_final=round_up_float32(max(raw_lower,0.72,raw_negative_upper*1.1,raw_negative_upper+0.05));all_live_finite;no_cap;accept_strictly_greater\n"+
			"blosum_final=round_up_float32(max(raw_lower,0,raw_negative_upper*1.1,raw_negative_upper+0.005));finite_only;accept_strictly_greater\n"+
			"retained_counts=original_calibration_guaranteed_minima_not_remeasured\n"+
			"retained_counts_scope=lower_cutoffs_only;final_cosine_and_blosum_retention_not_measured\n"+
			"retained_counts_precision=original_double_metrics;not_remeasured_for_float_query_metrics\n"+
			"query_scores=cascade_v2_float\nraw_nonfinite=reject\nbiological_recalibration=false\n");
		final int nonempty=rows.size()-zeros;
		write(out.resolve("summary.tsv"),SUMMARY_HEADER+"\n"+rows.size()+"\t"+nonempty+"\t"+zeros+"\t"+members+"\t"+
			nonempty+"\t"+disabled[0]+"\t"+disabled[1]+"\t"+disabled[2]+"\n");
		write(out.resolve("PASS"),PASS+"\n");
		final StringBuilder manifest=new StringBuilder("key\tvalue\nschema\t").append(FamilyGateCodec.SCHEMA)
			.append("\nraw_lower_sha80\t").append(lowerPin).append("\nraw_upper_sha80\t").append(upperPin).append('\n');
		for(String name : OUTPUT_FILES){manifest.append(name).append('\t').append(DigestSuffix.file(out.resolve(name).toString())).append('\n');}
		write(out.resolve("result.sha80.tsv"),manifest.toString());
		System.out.println(PASS+" families="+rows.size()+" nonempty="+nonempty+" zeros="+zeros+" members="+members);
	}

	static void readLower(Path path,ArrayList<Row> rows,HashMap<String,Row> byId,int limit){
		final ByteFile input=ByteFile.makeByteFile(path.toString(),false); final LineParser1 lp=new LineParser1('\t');
		try{
			header(input,FamilyCalibration.CUTOFF_HEADER);
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				lp.set(line); require(lp.terms()==10,"Wrong raw lower field count");
				final int rank=nonnegativeInt(lp,0), members=nonnegativeInt(lp,2); final String id=lp.parseString(1);
				require(rank>=0 && members>=0 && !id.isEmpty(),"Invalid raw lower identity/count");
				require(rows.size()<limit,"Raw lower exceeds maxfamilies");
				final Row row=new Row(rank,id,members);
				require(byId.put(id,row)==null,"Duplicate raw lower family: "+id);
				if(members==0){
					require(lp.termEquals("NO_MEMBERS",3),"Zero member status must be NO_MEMBERS");
					for(int field=4; field<=8; field+=2){
						require(lp.termEquals("NA",field) && nonnegativeInt(lp,field+1)==0,"Zero member lower must retain NA/0");
					}
				}else{
					require(lp.termEquals("CALIBRATED",3),"Live lower status must be CALIBRATED");
					for(int i=0; i<3; i++){
						row.lower[i]=Double.parseDouble(lp.parseString(4+2*i));
						domain(row.lower[i],i);
						row.retained[i]=nonnegativeInt(lp,5+2*i);
						require(row.retained[i]>=members-members/100 && row.retained[i]<=members,"Invalid raw retained count: "+id);
					}
				}
				rows.add(row);
			}
		}finally{require(!input.close(),"Raw lower read failed");}
		require(!rows.isEmpty(),"Raw lower must contain at least one family, including zero-member families");
	}

	static void readUpper(Path path,HashMap<String,Row> byId){readUpper(path, byId, 100000, 500);}
	static void readUpper(Path path,HashMap<String,Row> byId,int expectedNegativeCount,int expectedAlignedCount){
		require(expectedNegativeCount>0 && expectedAlignedCount>0 && expectedAlignedCount<=expectedNegativeCount,
			"Upper population must have positive aligned/query counts");
		final ByteFile input=ByteFile.makeByteFile(path.toString(),false); final LineParser1 lp=new LineParser1('\t');
		try{
			header(input,RAW_UPPER_HEADER);
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				lp.set(line); require(lp.terms()==7,"Wrong raw upper field count");
				final String id=lp.parseString(1); final Row row=byId.get(id);
				require(row!=null && row.rank==nonnegativeInt(lp,0) && row.members>0,"Unknown, mismatched or zero-member upper family: "+id);
				require(!row.hasUpper,"Duplicate raw upper family: "+id);
				row.negativeCount=nonnegativeInt(lp,2); row.alignedCount=nonnegativeInt(lp,4);
				require(row.negativeCount==expectedNegativeCount && row.alignedCount==expectedAlignedCount,
					"Unexpected raw upper population counts: "+id+" expected="+expectedNegativeCount+"/"+expectedAlignedCount+
					" observed="+row.negativeCount+"/"+row.alignedCount);
				for(int i=0; i<3; i++){row.upper[i]=Double.parseDouble(lp.parseString(i==0 ? 3 : i+4)); domain(row.upper[i],i);}
				row.hasUpper=true;
			}
		}finally{require(!input.close(),"Raw upper read failed");}
	}

	private static int nonnegativeInt(LineParser1 lp,int term){
		final int length=lp.setBounds(term);
		require(length>0 && length<=10,"Count/rank must fit a nonnegative int");
		for(int i=lp.a(); i<lp.b(); i++){require(lp.line()[i]>='0' && lp.line()[i]<='9',"Nondecimal count/rank");}
		// Ten digits fit long, whereas Parse.parseInt deliberately does not detect overflow.
		final long value=lp.parseLong(term);
		require(value<=Integer.MAX_VALUE,"Count/rank exceeds signed32 range");
		return (int)value;
	}

	private static void domain(double value,int metric){
		require(Double.isFinite(value),"Raw collection cutoff must be finite");
		require(metric==2 || (value>=0 && value<=(metric==0 ? 1 : 100)),"Raw cutoff outside metric domain");
	}

	private static void writeTables(Path out,ArrayList<Row> rows){
		final ByteStreamWriter low=writer(out.resolve("cutoffs.tsv")), high=writer(out.resolve("uppers.tsv"));
		final ByteBuilder a=new ByteBuilder(), b=new ByteBuilder();
		try{
			low.println(FamilyGateCodec.LOWER_HEADER); high.println(FamilyGateCodec.UPPER_HEADER);
			for(Row row : rows){
				a.clear().append(row.rank).tab().append(row.id).tab().append(row.members).tab();
				if(row.members==0){a.append("NO_MEMBERS\tNA\t0\tNA\t0\tNA\t0");}
				else{
					a.append("CALIBRATED");
					for(int i=0; i<3; i++){FamilyGateCodec.appendLower(a.tab(),row.lowerEncoded[i]).tab().append(row.retained[i]);}
					b.clear().append(row.rank).tab().append(row.id).tab().append(row.negativeCount).tab();
					FamilyGateCodec.appendCosineAccept(b,row.upperEncoded[0]).tab().append(row.alignedCount).tab();
					FamilyGateCodec.appendUpper(b,row.upperEncoded[1]).tab(); FamilyGateCodec.appendFinalBlosum(b,row.upperEncoded[2]);
					high.println(b);
				}
				low.println(a);
			}
		}finally{require(!(low.poisonAndWait() | high.poisonAndWait()),"Gate table writer failed");}
	}

	private static void pin(Path path,String expected){
		require(path.isAbsolute() && Files.isRegularFile(path) && !Files.isSymbolicLink(path),"Require regular absolute raw input: "+path);
		require(DigestSuffix.file(path.toString()).equals(expected),"Raw input pin mismatch: "+path);
	}
	private static void header(ByteFile file,String expected){require(Arrays.equals(file.nextLine(),expected.getBytes(StandardCharsets.US_ASCII)),"Unexpected raw gate header");}
	private static String required(HashMap<String,String> options,String key){final String v=options.get(key); require(v!=null && !v.isEmpty(),"Missing option: "+key); return v;}
	private static ByteStreamWriter writer(Path path){final ByteStreamWriter w=new ByteStreamWriter(path.toString(),false,false,false); w.start(); return w;}
	private static void write(Path path,String value) throws Exception{Files.write(path,value.getBytes(StandardCharsets.US_ASCII),java.nio.file.StandardOpenOption.CREATE_NEW);}
	private static void require(boolean condition,String message){if(!condition){throw new IllegalArgumentException(message);}}

	static final class Row{
		Row(int rank_,String id_,int members_){rank=rank_; id=id_; members=members_;}
		final int rank,members; final String id;
		final double[] lower=new double[3], upper=new double[3];
		final float[] lowerEncoded=new float[3], upperEncoded=new float[3];
		final int[] retained=new int[3];
		int negativeCount,alignedCount; boolean hasUpper;
	}
	public static final String RAW_UPPER_HEADER="roster_rank\tfamily_id\tnegative_count\tcosine_upper\taligned_count\tidentity_pct_upper\tblosum_a_sq_upper";
	public static final String SUMMARY_HEADER="families\tnonempty\tzeros\tmembers\tupper_rows\tdisabled_cosine\tdisabled_identity\tdisabled_blosum";
	public static final String PASS="FAMILY_GATE_EXPORT_PASS";
	private static final String[] OUTPUT_FILES={"cutoffs.tsv","uppers.tsv","recipe.txt","summary.tsv","PASS"};
}
