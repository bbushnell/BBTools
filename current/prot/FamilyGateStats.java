package prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/** Read-only threshold audit of pinned raw calibrations and a current uncapped gate export.
 * No sequence loading, recalibration, query scoring or gate rewriting.
 * @author Sayu
 */
public final class FamilyGateStats{

	private FamilyGateStats(){}
	public static void main(String[] args){
		try{
			final PreParser pp=new PreParser(args,FamilyGateStats.class,false);
			try{run(pp.args);}finally{Shared.closeStream(pp.outstream);}
		}catch(Throwable failure){failure.printStackTrace(); System.exit(1);}
	}

	public static void run(String[] args) throws Exception{
		final HashMap<String,String> options=new HashMap<String,String>();
		for(String arg : args){
			final String[] pair=Parse.splitOnFirst(arg,'=');
			require(pair.length==2,"Require flag=value arguments");
			final String key=pair[0].toLowerCase(Locale.ROOT);
			require(Arrays.asList("lowerindex","lowerindexsha80","upperindex","upperindexsha80",
				"gateindex","gateindexsha80","out","maxfamilies").contains(key) &&
				options.put(key,pair[1])==null,"Unknown/duplicate flag=value option");
		}
		final int limit=options.containsKey("maxfamilies") ? Parse.parseIntKMG(options.get("maxfamilies")) : 1000000;
		require(limit>0 && limit<=1000000,"maxfamilies must be in1..1000000");
		final Path out=absolute(required(options,"out"));
		require(!Files.exists(out) && !Files.isSymbolicLink(out) && out.getParent()!=null &&
			Files.isDirectory(out.getParent()) && out.getParent().toRealPath().equals(out.getParent()),"Require fresh output directory with real parent");
		final Inputs pins=new Inputs();
		final Bucket[][] indices=new Bucket[3][];
		for(int kind=0;kind<3;kind++){
			final String key=INDEX_NAMES[kind];
			final Path path=absolute(required(options,key)); pins.add(path,required(options,key+"sha80"));
			indices[kind]=index(path,kind);
		}
		final Row[] rows=new Row[limit]; final HashMap<String,Row> ids=new HashMap<String,Row>();
		int families=0,live=0,zeros=0,buckets=0; long members=0,residues=0;
		// Raw tables were produced/consumed with Java's correctly rounded decimal parser.
		// Use LineParser1's explicit Java fallback for bit parity, including exponent notation.
		final boolean oldParse=Tools.FORCE_JAVA_PARSE_DOUBLE; Tools.FORCE_JAVA_PARSE_DOUBLE=true;
		try{
			for(int b=0;b<256;b++){
				final Bucket lower=indices[0][b],upper=indices[1][b],gate=indices[2][b];
				if(lower==null){require(upper==null && gate==null,"Unexpected upper/gate bucket"); continue;}
				require(lower.same(upper) && lower.same(gate),"Index census mismatch at bucket "+b);
				buckets++; families=Math.addExact(families,lower.families); live=Math.addExact(live,lower.live);
				zeros=Math.addExact(zeros,lower.zeros); members=Math.addExact(members,lower.members); residues=Math.addExact(residues,lower.residues);
				require(families<=limit,"Family limit exceeded");
				pins.add(lower.manifest,lower.pin); pins.add(upper.manifest,upper.pin); pins.add(gate.manifest,gate.pin);
				final Path rawLower=rawPayload(lower.manifest,0,pins),rawUpper=rawPayload(upper.manifest,1,pins);
				final HashMap<String,String> gm=manifest(gate.manifest,"key\tvalue",GATE_KEYS);
				require(gm.size()==GATE_KEYS.length && FamilyGateCodec.SCHEMA.equals(gm.get("schema")),"Require complete current-formula final-gate manifest");
				require(pins.values.get(rawLower).equals(gm.get("raw_lower_sha80")) &&
					pins.values.get(rawUpper).equals(gm.get("raw_upper_sha80")),"Gate/raw provenance mismatch");
				for(String name : GATE_FILES){pins.add(gate.manifest.getParent().resolve(name),gm.get(name));}
				pass(gate.manifest.getParent().resolve("PASS"),FamilyGateExport.PASS);
				readTable(rawLower,0,lower,rows,ids);
				readTable(rawUpper,1,lower,rows,ids);
				readTable(gate.manifest.getParent().resolve("cutoffs.tsv"),2,lower,rows,ids);
				readTable(gate.manifest.getParent().resolve("uppers.tsv"),3,lower,rows,ids);
			}
		}finally{Tools.FORCE_JAVA_PARSE_DOUBLE=oldParse;}
		require(buckets>0 && families>0 && families==live+zeros && ids.size()==families,"Incomplete family census");
		for(int rank=0;rank<rows.length;rank++){
			final Row row=rows[rank];
			if(rank>=families){require(row==null,"Roster ranks must be contiguous"); continue;}
			require(row!=null && row.seen==(row.members==0 ? 5 : 15),"Missing joined family or gate row at rank "+rank);
		}
		pins.verify(); // All failures above occur before publishing any report.
		Files.createDirectory(out);
		writeReport(out,rows,families,live,zeros,buckets,members,residues);
		final StringBuilder inputText=new StringBuilder("path\tsha80\n");
		for(java.util.Map.Entry<Path,String> entry : pins.values.entrySet()){inputText.append(entry.getKey()).append('\t').append(entry.getValue()).append('\n');}
		write(out.resolve("inputs.tsv"),inputText.toString());
		write(out.resolve("recipe.txt"),"schema=family_gate_stats_v3\nformula=max(raw_lower,0,raw_upper*1.1,raw_upper+0.005)\n"+
			"quantile=nearest_rank;index=max(0,ceil(p*n)-1);no_interpolation\nempty_families=excluded_from_numeric_distributions\n"+
			"gap=raw_final-raw_lower;bins=0,(0,1e-6],(1e-6,1e-4],(1e-4,1e-3],(1e-3,1e-2],(1e-2,infinity)\n"+
			"determinants=LOWER:L>M;LOWER_TIE:L==M;otherwise_MARGIN\n"+
			"M=max(0,raw_upper*1.1,raw_upper+0.005);no_cap;0.144_counts_are_historical_reference_only\nrecalibration=false\nquery_replay=false\n");
		pins.verify();
		write(out.resolve("PASS"),PASS+"\n");
		final StringBuilder manifest=new StringBuilder("file\tsha80\n");
		for(String name : OUTPUTS){manifest.append(name).append('\t').append(DigestSuffix.file(out.resolve(name).toString())).append('\n');}
		write(out.resolve("result.sha80.tsv"),manifest.toString());
		System.out.println(PASS+" families="+families+" live="+live+" empty="+zeros);
	}

	private static Bucket[] index(Path path,int kind){
		final Bucket[] rows=new Bucket[256]; final ByteFile input=ByteFile.makeByteFile(path.toString(),false);
		final LineParser1 lp=new LineParser1('\t');
		try{
			header(input,INDEX_HEADERS[kind]);
			for(byte[] line=input.nextLine();line!=null;line=input.nextLine()){
				lp.set(line); require(lp.terms()==11-kind,"Wrong index field count");
				final int bucket=uint(lp,0); require(bucket<256 && rows[bucket]==null,"Duplicate/out-of-range bucket");
				final int f=uint(lp,1),n=uint(lp,2),z=uint(lp,3); final long m=ulong(lp,4),aa=ulong(lp,5);
				require((long)n+z==f && f<=1000000,"Invalid index census");
				if(kind==1){require(uint(lp,6)==n,"Upper population differs from live census");}
				final int column=8-kind;
				rows[bucket]=new Bucket(bucket,f,n,z,m,aa,absolute(lp.parseString(column)),DigestSuffix.requireSuffix(lp.parseString(column+1),"manifest"));
				require(lp.termEquals("NEW",lp.terms()-1) || lp.termEquals("REUSED_ACCEPTED",lp.terms()-1),"Unknown index lineage");
			}
		}finally{require(!input.close(),"Index read failed");}
		return rows;
	}

	private static HashMap<String,String> manifest(Path path,String expected,String[] keys){
		final HashMap<String,String> map=new HashMap<String,String>();
		final ByteFile input=ByteFile.makeByteFile(path.toString(),false); final LineParser1 lp=new LineParser1('\t');
		try{
			header(input,expected);
			for(byte[] line=input.nextLine();line!=null;line=input.nextLine()){
				lp.set(line); require(lp.terms()==2,"Malformed manifest");
				final String key=lp.parseString(0);
				require(Arrays.asList(keys).contains(key) && map.put(key,lp.parseString(1))==null,"Unknown/duplicate manifest key");
			}
		}finally{require(!input.close(),"Manifest read failed");}
		return map;
	}

	private static Path rawPayload(Path manifest,int kind,Inputs pins) throws Exception{
		final String[] names=kind==0 ? LOWER_FILES : UPPER_FILES;
		final HashMap<String,String> map=manifest(manifest,"file\tsha80",names);
		final String payload=kind==0 ? "cutoffs.tsv" : "uppers.tsv";
		final String chosen=map.containsKey(payload) ? payload : payload+".gz";
		final int expected=kind==0 ? 6 : 5;
		require(map.size()==expected && map.containsKey(chosen),"Incomplete raw manifest");
		for(String logical : kind==0 ? new String[]{"cutoffs.tsv","spectra.tsv","summary.tsv","recipe.txt","timing.tsv","PASS"} : UPPER_FILES){
			require(map.containsKey(logical) != map.containsKey(logical+".gz"),"Missing/duplicate logical raw artifact");
		}
		for(String hash : map.values()){DigestSuffix.requireSuffix(hash,"raw artifact");}
		final Path path=manifest.getParent().resolve(chosen); pins.add(path,map.get(chosen));
		final Path pass=manifest.getParent().resolve("PASS"); pins.add(pass,map.get("PASS"));
		pass(pass,kind==0 ? "FAMILY_ALL_MEMBER_CALIBRATION_PASS" : "FAMILY_NEGATIVE_UPPER_CALIBRATION_PASS");
		return path; // Unconsumed spectra/other bulk payloads are intentionally not opened.
	}

	private static void readTable(Path path,int kind,Bucket bucket,Row[] rows,HashMap<String,Row> ids){
		final ByteFile input=ByteFile.makeByteFile(path.toString(),false); final LineParser1 lp=new LineParser1('\t');
		int count=0,live=0; long members=0;
		try{
			header(input,TABLE_HEADERS[kind]);
			for(byte[] line=input.nextLine();line!=null;line=input.nextLine()){
				lp.set(line); require(lp.terms()==(kind%2==0 ? 10 : 7),"Wrong cutoff table field count");
				final int rank=uint(lp,0); require(rank<rows.length,"Rank exceeds bound");
				final String id=lp.parseString(1); require(!id.isEmpty(),"Empty family ID");
				Row row=rows[rank];
				if(kind==0){
					require(row==null,"Duplicate raw lower rank"); row=new Row(rank,id,uint(lp,2),bucket.bucket); rows[rank]=row;
					require(ids.put(id,row)==null,"Duplicate raw lower family ID");
				}else{require(row!=null && row.id.equals(id) && row.bucket==bucket.bucket,"Unknown/cross-bucket/mismatched family");}
				require((row.seen&(1<<kind))==0,"Duplicate cutoff row"); row.seen|=1<<kind;
				if(kind%2==0){
					require(uint(lp,2)==row.members,"Lower member mismatch"); members=Math.addExact(members,row.members);
					if(row.members==0){
						require(lp.termEquals("NO_MEMBERS",3),"Invalid empty status");
						for(int i=4;i<=8;i+=2){require(lp.termEquals("NA",i) && uint(lp,i+1)==0,"Empty lower must be NA/0");}
					}else{
						live++; require(lp.termEquals("CALIBRATED",3),"Invalid live status");
						for(int i=0;i<3;i++){
							final int retained=uint(lp,5+2*i);
							require(retained>=row.members-row.members/100 && retained<=row.members,"Invalid retained minimum");
							if(kind==0){row.lower[i]=decimal(lp,4+2*i,i); row.retained[i]=retained;}
							else{require(retained==row.retained[i] && bits(FamilyGateCodec.readLower(lp,4+2*i))==bits(FamilyGateCodec.lower(row.lower[i])),"Encoded lower differs from raw");}
						}
					}
				}else{
					require(row.members>0 && uint(lp,2)==100000 && uint(lp,4)==500,"Unexpected upper population");
					if(kind==1){for(int i=0;i<3;i++){row.upper[i]=decimal(lp,i==0 ? 3 : i+4,i);}}
					else{
						final float cosine=FamilyGateCodec.readCosineAccept(lp,3);
						final float expectedCosine=roundUp(Math.max(row.lower[0],Math.max(0.72,Math.max(row.upper[0]*1.1,row.upper[0]+0.05))));
						require(bits(cosine)==bits(expectedCosine) && bits(cosine)==bits(FamilyGateCodec.cosineAccept(row.lower[0],row.upper[0])),"Encoded cosine differs from original-double margin");
						require(bits(FamilyGateCodec.readUpper(lp,5))==bits(FamilyGateCodec.upper(row.lower[1],row.upper[1])),"Encoded identity shortcut differs from raw");
						row.encoded=FamilyGateCodec.readFinalBlosum(lp,6);
						final double raw=rawFinal(row.lower[2],row.upper[2]); final float expected=roundUp(raw);
						require(bits(expected)==bits(FamilyGateCodec.finalBlosum(row.lower[2],row.upper[2])) && bits(expected)==bits(row.encoded),"Final gate differs from independent formula/rounding");
					}
				}
				count++;
			}
		}finally{require(!input.close(),"Cutoff table read failed");}
		require(count==(kind%2==0 ? bucket.families : bucket.live),"Cutoff row census mismatch");
		if(kind%2==0){require(live==bucket.live && members==bucket.members,"Lower live/member census mismatch");}
	}

	/** Independent formula oracle: no codec call, no rounded lower input. */
	public static double rawFinal(double lower,double upper){
		require(Double.isFinite(lower) && Double.isFinite(upper) && Math.abs(lower)<=Float.MAX_VALUE && Math.abs(upper)<=Float.MAX_VALUE,"Raw BLOSUM must fit finite float32");
		return Math.max(lower,margin(upper)); // Math.max preserves the specified signed-zero behavior.
	}
	public static float roundUp(double value){
		require(Double.isFinite(value) && Math.abs(value)<=Float.MAX_VALUE,"Nonrepresentable final threshold");
		final float nearest=(float)value; return nearest<value ? Math.nextUp(nearest) : nearest;
	}
	public static String determinant(double lower,double upper){
		rawFinal(lower,upper); final double margin=margin(upper);
		return lower>margin ? "LOWER" : lower==margin ? "LOWER_TIE" : "MARGIN";
	}
	private static double margin(double upper){return Math.max(0,Math.max(upper*1.1,upper+0.005));}
	/** Sorted finite input, exact integer nearest-rank index; empty summaries are handled by the caller. */
	public static double quantile(double[] sorted,int percent){
		require(sorted.length>0 && percent>=0 && percent<=100,"Invalid quantile request");
		final int index=(int)Math.max(0,((long)percent*sorted.length+99)/100-1);
		assert(index<sorted.length) : "Nearest-rank percentile must select one observed threshold, not interpolate past the sample";
		return sorted[index];
	}
	public static int gapBin(double gap){
		require(Double.isFinite(gap) && gap>=0,"Final raw threshold cannot be below raw lower");
		if(gap==0){return 0;}
		for(int i=0;i<GAP_LIMITS.length;i++){if(gap<=GAP_LIMITS[i]){return i+1;}}
		return 5;
	}

	private static void writeReport(Path out,Row[] rows,int families,int live,int zeros,int buckets,long members,long residues) throws Exception{
		final double[][] data=new double[7][live]; final long[] categories=new long[CATEGORIES.length],gaps=new long[6];
		long above=0,rawCap=0,encodedDecimalCap=0,encodedCapBits=0; int slot=0;
		final ByteStreamWriter writer=new ByteStreamWriter(out.resolve("families.tsv").toString(),false,false,false); writer.start();
		final ByteBuilder bb=new ByteBuilder();
		final double[] values=new double[7];
		try{
			writer.println("roster_rank\tfamily_id\tmembers\tstatus\traw_lower\traw_upper\tmargin\tpolicy_margin\traw_final\tencoded_final\tencoded_bits_u32\traw_gap\tencoded_gap\trounding_gap\tdeterminant");
			for(int rank=0;rank<families;rank++){
				final Row row=rows[rank]; bb.clear().append(rank).tab().append(row.id).tab().append(row.members).tab();
				if(row.members==0){writer.println(bb.append("NO_MEMBERS\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNO_MEMBERS")); continue;}
				final double l=row.lower[2],u=row.upper[2],f=rawFinal(l,u),e=row.encoded;
				final String category=determinant(l,u); categories[Arrays.asList(CATEGORIES).indexOf(category)]++; gaps[gapBin(f-l)]++;
				if(l>0.144){above++;} if(f==0.144){rawCap++;} if(e==0.144){encodedDecimalCap++;}
				if(bits(row.encoded)==bits(roundUp(0.144))){encodedCapBits++;}
				values[0]=l; values[1]=u; values[2]=f; values[3]=e; values[4]=f-l; values[5]=e-l; values[6]=e-f;
				for(int i=0;i<values.length;i++){require(Double.isFinite(values[i]),"Nonfinite report statistic"); data[i][slot]=values[i];} slot++;
				bb.append("CALIBRATED");
				bb.tab().append(Double.toString(l)).tab().append(Double.toString(u)).tab().append(Double.toString(margin(u)))
					.tab().append(Double.toString(margin(u))).tab().append(Double.toString(f)).tab().append(Double.toString(e));
				bb.tab().append(Integer.toUnsignedLong(bits(row.encoded)));
				for(int i=4;i<7;i++){bb.tab().append(Double.toString(values[i]));}
				writer.println(bb.tab().append(category));
			}
		}finally{require(!writer.poisonAndWait(),"Report writer failed");}
		require(slot==live,"Numeric count must exclude exactly the empty families");
		final StringBuilder quantiles=new StringBuilder("metric\tn\tmin\tp10\tp25\tmedian\tp75\tp90\tmax\n");
		for(int i=0;i<data.length;i++){
			Arrays.sort(data[i]); quantiles.append(METRICS[i]).append('\t').append(live);
			for(int p : new int[]{0,10,25,50,75,90,100}){quantiles.append('\t').append(live==0 ? "NA" : Double.toString(quantile(data[i],p)));}
			quantiles.append('\n');
		}
		write(out.resolve("quantiles.tsv"),quantiles.toString());
		final StringBuilder summary=new StringBuilder("metric\tcount\n");
		final String[] names={"buckets","families","live","empty","members","residues","raw_lower_above_0144","raw_final_equal_0144","encoded_equal_decimal_0144","encoded_equal_upward_0144_bits"};
		final long[] counts={buckets,families,live,zeros,members,residues,above,rawCap,encodedDecimalCap,encodedCapBits};
		for(int i=0;i<names.length;i++){summary.append(names[i]).append('\t').append(counts[i]).append('\n');}
		long classified=0,binned=0;
		for(long count : categories){classified+=count;} for(long count : gaps){binned+=count;}
		require(classified==live && binned==live,"Every live family needs exactly one determinant and gap bin");
		for(int i=0;i<CATEGORIES.length;i++){summary.append(CATEGORIES[i]).append('\t').append(categories[i]).append('\n');}
		write(out.resolve("summary.tsv"),summary.toString());
		final StringBuilder histogram=new StringBuilder("raw_final_minus_lower_bin\tcount\n");
		for(int i=0;i<gaps.length;i++){histogram.append(GAP_NAMES[i]).append('\t').append(gaps[i]).append('\n');}
		write(out.resolve("gaps.tsv"),histogram.toString());
	}

	private static double decimal(LineParser1 lp,int term,int metric){
		final double value=lp.parseDouble(term);
		require(Double.isFinite(value) && Math.abs(value)<=Float.MAX_VALUE && (metric==2 || value>=0 && value<=(metric==0 ? 1 : 100)),"Invalid raw cutoff");
		return value;
	}
	private static long ulong(LineParser1 lp,int term){
		lp.setBounds(term); require(lp.b()>lp.a(),"Empty integer"); long value=0;
		for(int i=lp.a();i<lp.b();i++){final int digit=lp.line()[i]-'0'; require(digit>=0 && digit<=9 && value<=(Long.MAX_VALUE-digit)/10,"Nondecimal/overflowed integer"); value=value*10+digit;}
		return value;
	}
	private static int uint(LineParser1 lp,int term){final long value=ulong(lp,term); require(value<=Integer.MAX_VALUE,"Integer exceeds signed32"); return (int)value;}
	private static int bits(float f){return Float.floatToRawIntBits(f);}
	private static Path absolute(String value){final Path p=Paths.get(value); require(p.isAbsolute() && p.normalize().equals(p),"Require normalized absolute path"); return p;}
	private static void pass(Path p,String text) throws Exception{require(Arrays.equals(Files.readAllBytes(p),(text+"\n").getBytes(StandardCharsets.US_ASCII)),"Missing producer PASS");}
	private static void header(ByteFile input,String text){require(Arrays.equals(input.nextLine(),text.getBytes(StandardCharsets.US_ASCII)),"Unexpected header: "+text);}
	private static void write(Path p,String text) throws Exception{Files.write(p,text.getBytes(StandardCharsets.UTF_8),StandardOpenOption.CREATE_NEW);}
	private static String required(HashMap<String,String> map,String key){final String v=map.get(key); require(v!=null && !v.isEmpty(),"Missing "+key); return v;}
	private static void require(boolean condition,String why){if(!condition){throw new IllegalArgumentException(why);}}

	private static final class Inputs{
		void add(Path p,String hash) throws Exception{
			DigestSuffix.requireSuffix(hash,"input"); require(values.size()<4096,"Too many report inputs");
			final String old=values.put(p,hash); require(old==null || old.equals(hash),"Conflicting pins for one path"); check(p,hash);
		}
		void verify() throws Exception{for(java.util.Map.Entry<Path,String> e : values.entrySet()){check(e.getKey(),e.getValue());}}
		void check(Path p,String hash) throws Exception{require(Files.isRegularFile(p) && !Files.isSymbolicLink(p) && p.toRealPath().equals(p) && DigestSuffix.file(p.toString()).equals(hash),"Input path/pin mismatch: "+p);}
		final LinkedHashMap<Path,String> values=new LinkedHashMap<Path,String>();
	}
	private static final class Bucket{
		Bucket(int b,int f,int n,int z,long m,long aa,Path p,String h){bucket=b; families=f; live=n; zeros=z; members=m; residues=aa; manifest=p; pin=h;}
		boolean same(Bucket b){return b!=null && bucket==b.bucket && families==b.families && live==b.live && zeros==b.zeros && members==b.members && residues==b.residues;}
		final int bucket,families,live,zeros; final long members,residues; final Path manifest; final String pin;
	}
	private static final class Row{
		Row(int r,String i,int m,int b){rank=r; id=i; members=m; bucket=b;}
		final int rank,members,bucket; final String id; final double[] lower=new double[3],upper=new double[3]; final int[] retained=new int[3]; int seen; float encoded;
	}
	public static final String PASS="FAMILY_GATE_STATS_PASS";
	public static final String[] INDEX_HEADERS={
		"bucket\tfamilies\tnonempty\tzeros\tmembers\tresidues\treferences_manifest\treferences_sha80\tlower_manifest\tlower_sha80\tlineage",
		"bucket\tfamilies\tnonempty\tzeros\tmembers\tresidues\tupper_rows\tupper_manifest\tupper_sha80\tlineage",
		"bucket\tfamilies\tnonempty\tzeros\tmembers\tresidues\tgate_manifest\tgate_sha80\tlineage"};
	private static final String[] INDEX_NAMES={"lowerindex","upperindex","gateindex"};
	private static final String[] TABLE_HEADERS={FamilyCalibration.CUTOFF_HEADER,FamilyGateExport.RAW_UPPER_HEADER,FamilyGateCodec.LOWER_HEADER,FamilyGateCodec.UPPER_HEADER};
	private static final String[] LOWER_FILES={"cutoffs.tsv","cutoffs.tsv.gz","spectra.tsv","spectra.tsv.gz","summary.tsv","recipe.txt","timing.tsv","PASS"};
	private static final String[] UPPER_FILES={"uppers.tsv","summary.tsv","timing.tsv","recipe.txt","PASS"};
	private static final String[] GATE_KEYS={"schema","raw_lower_sha80","raw_upper_sha80","cutoffs.tsv","uppers.tsv","recipe.txt","summary.tsv","PASS"};
	private static final String[] GATE_FILES={"cutoffs.tsv","uppers.tsv","recipe.txt","summary.tsv","PASS"};
	private static final String[] OUTPUTS={"families.tsv","quantiles.tsv","summary.tsv","gaps.tsv","inputs.tsv","recipe.txt","PASS"};
	private static final String[] METRICS={"raw_lower","raw_upper","raw_final","encoded_final","raw_gap","encoded_gap","rounding_gap"};
	private static final String[] CATEGORIES={"LOWER","LOWER_TIE","MARGIN"};
	private static final double[] GAP_LIMITS={1e-6,1e-4,1e-3,1e-2};
	private static final String[] GAP_NAMES={"0","(0,1e-6]","(1e-6,1e-4]","(1e-4,1e-3]","(1e-3,1e-2]","(1e-2,infinity)"};
}
