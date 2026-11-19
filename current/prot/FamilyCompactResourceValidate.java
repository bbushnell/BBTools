package prot;

import java.io.BufferedWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import ddl.DDLFamilyProfiles;
import ddl.DDLPentamerBlacklist;
import ddl.ProteinKmerGeometry;
import fileIO.ByteFile;
import parse.LineParser1;
import parse.Parse;

/** Production compact-resource equality validator and small binding helper.
 * @author Yelan
 */
public final class FamilyCompactResourceValidate{

	private FamilyCompactResourceValidate(){}

	public static void main(String[] args) throws Exception{
		final HashMap<String,String> options=parse(args);
		final String mode=options.containsKey("mode") ? options.get("mode") : "compare";
		if("matrixbinding".equals(mode)){
			final String path=required(options,"matrixresource");
			final ProteinSubstitutionMatrix matrix=ProteinSubstitutionMatrix.load(
				ProteinSubstitutionMatrix.POSSUM55_NAME,path);
			if(options.containsKey("acceptedmatrix")){
				final ProteinSubstitutionMatrix accepted=ProteinSubstitutionMatrix.load(
					"custom40",required(options,"acceptedmatrix"));
				require(scoreBinding(matrix).equals(scoreBinding(accepted)),
					"POSSUM55 projection must alter only matrix name/description");
			}
			System.out.println(matrix.binding());
			return;
		}
		if("blacklist".equals(mode)){
			final Path out=path(required(options,"out"));
			require(!Files.exists(out) && !Files.isSymbolicLink(out) && out.getParent()!=null &&
				Files.isDirectory(out.getParent()),"Require a fresh output file with an existing parent");
			final BlacklistFasta shortList=readBlacklist(path(required(options,"short")));
			final BlacklistFasta tailsList=readBlacklist(path(required(options,"tails")));
			final Counts counts=writeNeutralBlacklist(shortList,tailsList,out);
			final String pass="FAMILY_COMPACT_BLACKLIST_PASS words="+counts.physical+
				" k="+shortList.header.k+" mode="+shortList.header.mode+"\n";
			System.out.print(pass);
			return;
		}
		require("compare".equals(mode),"Unknown mode: "+mode);
		final Path out=path(required(options,"out"));
		require(!Files.exists(out) && !Files.isSymbolicLink(out) && out.getParent()!=null &&
			Files.isDirectory(out.getParent()),"Require a fresh output directory with an existing parent");
		final String view=required(options,"view");
		final Expected expected=expected(view,options);
		final int threads=integer(options,"threads",1),maxFamilies=integer(options,"maxfamilies",1000000);
		require(threads>0 && maxFamilies>0,"threads and maxfamilies must be positive");
		final FamilyReferenceLoader.Loaded legacy=FamilyReferenceLoader.loadMetricBound(
			required(options,"references"),threads,maxFamilies);
		final FamilyCompactLoader.Loaded compact=FamilyCompactLoader.load(
			required(options,"compact"),threads,maxFamilies);
		final char[][] profiles=readProfiles(required(options,"profilevectors"),legacy.all);
		final Sidecars sidecars=sidecars(legacy.all,required(options,"lengthmedians"),
			required(options,"familyhydro"),required(options,"familygc"),required(options,"familygcmap"),expected);
		final Counts counts=compare(expected,legacy,compact,profiles,sidecars);
		Files.createDirectory(out);
		Files.write(out.resolve("equality.tsv"),counts.tsv().getBytes(StandardCharsets.US_ASCII));
		final String pass="FAMILY_COMPACT_RESOURCE_VALIDATE_PASS view="+view+
			" physical="+counts.physical+" live="+counts.live+" members="+counts.members+
			" residues="+counts.residues+" valid_words="+counts.validWords+"\n";
		Files.write(out.resolve("PASS"),pass.getBytes(StandardCharsets.US_ASCII));
		System.out.print(pass);
	}

	private static Counts compare(Expected expected,FamilyReferenceLoader.Loaded legacy,FamilyCompactLoader.Loaded compact,
			char[][] profiles,Sidecars sidecars){
		validateHeader(expected,compact.header);
		require(legacy.all.length==compact.all.length &&
			legacy.live.length==compact.live.length,"Physical/live census mismatch");
		require(legacy.members==compact.members && legacy.residues==compact.residues &&
			legacy.validWords==compact.validWords,"Reference totals mismatch");
		long members=0,residues=0,validWords=0; int live=0;
		for(int i=0;i<legacy.all.length;i++){
			final FamilyReference a=legacy.all[i],b=compact.all[i];
			require(a.rank==i && b.rank==i,"References must be dense physical-rank ordered at ordinal "+i);
			equal(a.id,b.id,"family id",i);
			require(a.members==b.members && a.residues==b.residues && a.validWords==b.validWords,
				"Count mismatch at "+a.id);
			require(raw(a.cosineLower)==raw(b.cosineLower) && raw(a.identityLower)==raw(b.identityLower) &&
				raw(a.blosumLower)==raw(b.blosumLower) && raw(a.cosineAccept)==raw(b.cosineAccept) &&
				raw(a.identityUpper)==raw(b.identityUpper) && raw(a.blosumFinal)==raw(b.blosumFinal),
				"Gate raw32 mismatch at "+a.id);
			require(a.scalarMask==b.scalarMask && rawEquals(a.scalarBounds,b.scalarBounds),
				"Scalar raw64 mismatch at "+a.id);
			require(a.alignmentCalibration==b.alignmentCalibration,
				"Alignment calibration mismatch at "+a.id);
			if(a.members==0){
				require(compact.allStableIds[i]==i,"Empty stable ID mismatch at "+a.id);
			}else{
				require(Arrays.equals(a.consensusRaw,b.consensusRaw),"Consensus mismatch at "+a.id);
				require(Arrays.equals(profiles[i],b.ddl),"DDL mismatch at "+a.id);
				require(rawEquals(a.spectrum,b.spectrum),"Monomer mismatch at "+a.id);
				require(compact.live[live]==b && compact.liveStableIds[live]==sidecars.stableIds[i],
					"Live order or stable ID mismatch at "+a.id);
				require(raw(compact.lengthMedians[live])==raw(sidecars.lengthMedians[i]) &&
					raw(compact.kdMeans[live])==raw(sidecars.kdMeans[i]) &&
					compact.gcMeasured[live]==a.members &&
					raw(compact.gcMeans[live])==raw(sidecars.gcMeans[i]) &&
					raw(compact.gcSampleSds[live])==raw(sidecars.gcSampleSds[i]),
					"Sidecar raw64 mismatch at "+a.id);
				live++;
			}
			members=Math.addExact(members,a.members);
			residues=Math.addExact(residues,a.residues);
			validWords=Math.addExact(validWords,a.validWords);
		}
		require(live==legacy.live.length,"Live comparison count mismatch");
		require(legacy.all.length==expected.families && live==expected.live &&
			legacy.all.length-live==expected.empty && members==expected.members,
			"Accepted compact population constants mismatch");
		return new Counts(legacy.all.length,live,members,residues,validWords);
	}

	private static void validateHeader(Expected expected,FamilyCompactCodec.Header header){
		require(expected.view.equals(header.get("view")),"Compact header view mismatch");
		requireHeader(header,"families",Integer.toString(expected.families));
		requireHeader(header,"live",Integer.toString(expected.live));
		requireHeader(header,"empty",Integer.toString(expected.empty));
		require("alphabet.config".equals(header.get("alphabet")) &&
			"blacklist.faa.gz".equals(header.get("blacklist")) &&
			"POSSUM55.tsv".equals(header.get("matrixresource")),
			"Compact header shared-resource locators mismatch");
		require("512".equals(header.get("width")) && "6".equals(header.get("k")) &&
			"positional-v1".equals(header.get("sketch")) && "combined25".equals(header.get("rankingfeatures")),
			"Compact geometry or ranking header mismatch");
		requireHeader(header,"seed","0");
		requireHeader(header,"namespace","tail");
		requireHeader(header,"mode","nonempty");
		requireHeader(header,"exponent","6");
		requireHeader(header,"geometry",expected.geometry.name());
		requireHeader(header,"partitions","4");
		requireHeader(header,"overlap","2");
		requireHeader(header,"tag_bits","2");
		requireHeader(header,"segment_words","49");
		requireHeader(header,"formulaVersion",FamilyReferenceCodec.FORMULA_VERSION);
		requireHeader(header,"cosine","ranking");
		require(ProteinSubstitutionMatrix.POSSUM55_NAME.equals(header.get("matrix")) &&
			expected.matrixBinding.equals(header.get("matrixbinding")) &&
			"POSSUM55".equals(header.get("matrixdisplay")) &&
			ProteinSubstitutionMatrix.POSSUM55_DESCRIPTION.equals(header.get("matrixdescription")) &&
			"bbgap_possum55_20261116_v1".equals(header.get("matrixprovenance")) &&
			header.get("matrixbinding")!=null,
			"Compact POSSUM55 header mismatch");
		requireHeader(header,"top","50");
		requireHeader(header,"ddlreject","6");
		requireHeader(header,"ddlaccept","1025");
		requireHeader(header,"indexearlyreject","false");
		requireHeader(header,"scalargates","true");
		requireHeader(header,"scalarscale","0.1");
		requireHeader(header,"scalarz","3.8906");
		requireHeader(header,"identityacceptfloor","30");
		requireHeader(header,"blosumacceptfloor","0.13");
		for(Map.Entry<String,String> e:expected.scalarBounds.entrySet()){
			require(e.getValue().equals(header.get(e.getKey())),
				"Compact scalar bound mismatch: "+e.getKey());
		}
		require("reduced_literal".equals(header.get("blacklistmode")) &&
			"6".equals(header.get("blacklistk")) && "2048".equals(header.get("blacklistwords")) &&
			"bbgap_blacklist_k6_2048_20261116_v1".equals(header.get("blacklistprovenance")),
			"Compact blacklist header mismatch");
		requireHeader(header,"alignmentmode","pin_right_trim_left_aligned_query");
		requireHeader(header,"alignmentminspan",Integer.toString(expected.alignmentMinSpan));
		requireHeader(header,"alignmentminpaired",Integer.toString(expected.alignmentMinPaired));
		requireHeader(header,"heuristic","hitfourth_h020_c150");
		requireHeader(header,"minqueryaa",Integer.toString(expected.minAa));
		requireHeader(header,"maxqueryaa",Integer.toString(expected.maxAa));
		requireHeader(header,"minmemberaa",Integer.toString(expected.minAa));
		requireHeader(header,"maxmemberaa",Integer.toString(expected.maxAa));
		requireHeader(header,"taillow",Integer.toString(expected.minAa));
		requireHeader(header,"tailhigh",Integer.toString(expected.maxAa));
		requireHeader(header,"marginscale",expected.marginScale);
		requireHeader(header,"familygclambda","24");
		requireHeader(header,"gcfilter","true");
		requireHeader(header,"gc_z","8");
		requireHeader(header,"gc_minimum_members","30");
		requireHeader(header,"familygcfamilies",Integer.toString(expected.live));
		requireHeader(header,"familygcmembers",Long.toString(expected.members));
		requireHeader(header,"packageprovenance","bbgap_release_20261116_v1");
	}

	private static void requireHeader(FamilyCompactCodec.Header header,String key,String expected){
		require(expected.equals(header.get(key)),"Compact header "+key+" mismatch");
	}

	private static Sidecars sidecars(FamilyReference[] refs,String lengths,String hydro,
			String gc,String map,Expected expected) throws Exception{
		final HashMap<String,FamilyReference> byId=new HashMap<String,FamilyReference>(refs.length*2+1);
		for(int i=0;i<refs.length;i++){
			require(refs[i].rank==i,"References must be rank ordered before sidecar comparison");
			require(byId.put(refs[i].id,refs[i])==null,"Duplicate reference ID: "+refs[i].id);
		}
		final Sidecars s=new Sidecars(refs.length);
		readLengths(refs,byId,Paths.get(lengths),s);
		readHydro(refs,Paths.get(hydro),s);
		readGc(refs,Paths.get(gc),Paths.get(map),s,expected);
		return s;
	}

	private static void readLengths(FamilyReference[] refs,Map<String,FamilyReference> byId,
			Path path,Sidecars s) throws Exception{
		final ByteFile in=ByteFile.makeByteFile(path.toString(),false);
		final LineParser1 lp=new LineParser1('\t');
		int live=0;
		try{
			header(in,"family\tmembers\tlength_p10\tlength_median\tlength_p90\tcount_under_200");
			for(byte[] line=in.nextLine();line!=null;line=in.nextLine()){
				lp.set(line); require(lp.terms()==6,"Length median field count mismatch");
				final FamilyReference f=byId.get(lp.parseString(0));
				require(f!=null && lp.parseLong(1)==f.members,"Length median identity/count mismatch");
				if(f.members==0){require(lp.termEquals("NA",3),"Empty length median must be NA");}
				else{s.lengthMedians[f.rank]=slowDouble(lp,3); s.seenLength[f.rank]=true; live++;}
			}
		}finally{require(!in.close(),"Length median read failed");}
		require(live==live(refs,s.seenLength),"Length medians must cover every live family");
	}

	private static void readHydro(FamilyReference[] refs,Path path,Sidecars s) throws Exception{
		final ByteFile in=ByteFile.makeByteFile(path.toString(),false);
		final LineParser1 lp=new LineParser1('\t');
		int live=0,rows=0;
		try{
			header(in,"family_rank\tfamily_id\tN\thydrophobicity_mean");
			for(byte[] line=in.nextLine();line!=null;line=in.nextLine()){
				lp.set(line); require(lp.terms()==4,"Hydrophobicity field count mismatch");
				final int rank=lp.parseInt(0);
				require(rank==rows && rank<refs.length,"Hydrophobicity rank out of range: "+rank);
				final FamilyReference f=refs[rank];
				require(f.id.equals(lp.parseString(1)) && lp.parseLong(2)==f.members,
					"Hydrophobicity identity/count mismatch");
				if(f.members==0){require(lp.termEquals("NA",3),"Empty hydrophobicity must be NA");}
				else{s.kdMeans[rank]=slowDouble(lp,3); s.seenHydro[rank]=true; live++;}
				rows++;
			}
		}finally{require(!in.close(),"Hydrophobicity read failed");}
		require(rows==refs.length,"Hydrophobicity table must cover the physical roster");
		require(live==live(refs,s.seenHydro),"Hydrophobicities must cover every live family");
	}

	private static void readGc(FamilyReference[] refs,Path report,Path map,Sidecars s,
			Expected expected) throws Exception{
		final ByteFile reports=ByteFile.makeByteFile(report.toString(),false);
		final ByteFile maps=ByteFile.makeByteFile(map.toString(),false);
		final LineParser1 r=new LineParser1('\t'),m=new LineParser1('\t');
		try{
			header(reports,FamilyGcInputs.TAIL_HEADER_TEXT);
			header(maps,"stable_family_id\troster_rank\tfamily_id\tassigned\tresidues");
			for(byte[] rline=reports.nextLine(),mline=maps.nextLine();rline!=null || mline!=null;
					rline=reports.nextLine(),mline=maps.nextLine()){
				require(rline!=null && mline!=null,"GC report/map row count mismatch");
				r.set(rline); m.set(mline);
				require(r.terms()==10 && m.terms()==5,"GC field count mismatch");
				final int stable=r.parseInt(0),rank=r.parseInt(1);
				require(rank>=0 && rank<refs.length && m.parseInt(0)==stable &&
					m.parseInt(1)==rank,"GC rank/stable mismatch");
				final FamilyReference f=refs[rank];
				require(f.members>0 && !s.seenGc[rank] && f.id.equals(r.parseString(2)) &&
					f.id.equals(m.parseString(2)) && r.parseLong(3)==f.members &&
					r.parseLong(5)==f.members && m.parseLong(3)==f.members &&
					r.parseLong(4)==f.residues && m.parseLong(4)==f.residues,
					"GC identity/count mismatch");
				require(r.parseInt(8)==expected.minAa && r.parseInt(9)==expected.maxAa,
					"GC source tail interval mismatch");
				s.stableIds[rank]=stable;
				s.gcMeans[rank]=slowDouble(r,6);
				s.gcSampleSds[rank]=slowDouble(r,7);
				s.seenGc[rank]=true;
			}
		}finally{
			final boolean reportsFailed=reports.close(),mapsFailed=maps.close();
			require(!reportsFailed && !mapsFailed,"GC report/map read failed");
		}
		live(refs,s.seenGc);
	}

	private static int live(FamilyReference[] refs,boolean[] seen){
		int live=0;
		for(FamilyReference f:refs){
			if(f.members>0){require(seen[f.rank],"Missing live sidecar row for "+f.id);live++;}
		}
		return live;
	}

	/** Uses Java decimal rounding and exponent support to match exported raw binary64 values. */
	private static double slowDouble(LineParser1 lp,int field){
		lp.setBounds(field);
		return Parse.parseDoubleSlow(lp.line(),lp.a(),lp.b());
	}

	private static HashMap<String,String> parse(String[] args){
		final HashMap<String,String> map=new HashMap<String,String>();
		for(String arg:args){
			final int equals=arg.indexOf('=');
			require(equals>0 && equals<arg.length()-1,"Expected key=value argument: "+arg);
			String key=arg.substring(0,equals).toLowerCase();
			if("t".equals(key)){key="threads";}
			require(map.put(key,arg.substring(equals+1))==null,"Duplicate argument: "+key);
		}
		return map;
	}

	private static Path path(String text){return Paths.get(text).toAbsolutePath().normalize();}
	private static String required(HashMap<String,String> options,String key){
		final String value=options.get(key);
		require(value!=null && !value.isEmpty(),"Missing required argument: "+key);
		return value;
	}
	private static int integer(HashMap<String,String> options,String key,int defaultValue){
		final String value=options.get(key);
		return value==null ? defaultValue : Integer.parseInt(value);
	}
	private static Expected expected(String view,HashMap<String,String> options) throws Exception{
		final Expected e;
		if("short".equals(view)){e=new Expected(view,159994,159994,0,11149289,40,99,28,20,"0.985");}
		else if("tails".equals(view)){e=new Expected(view,120000,120000,0,117316888,100,100,60,43,"1");}
		else{
			require("fixture".equals(view),"Unknown compact validation view: "+view);
			e=new Expected(view,integer(options,"families",-1),integer(options,"live",-1),
				integer(options,"empty",-1),Long.parseLong(required(options,"members")),
				integer(options,"minaa",-1),integer(options,"maxaa",-1),
				integer(options,"alignmentminspan",-1),integer(options,"alignmentminpaired",-1),
				required(options,"marginscale"));
			require(e.families>=0 && e.live>=0 && e.empty>=0 && e.families==e.live+e.empty &&
				e.members>=0 && e.minAa>0 && e.maxAa>=e.minAa &&
				e.alignmentMinSpan>0 && e.alignmentMinPaired>0,
				"Custom compact validation views require explicit nonnegative census and positive range/span");
		}
		e.geometry=ProteinKmerGeometry.make(6,readAabits(path(required(options,"alphabet"))));
		e.matrixBinding=required(options,"matrixbinding");
		readScalarBounds(path(required(options,"scalarbounds")),e);
		return e;
	}
	private static String readAabits(Path path) throws Exception{
		final ByteFile in=ByteFile.makeByteFile(path.toString(),false);
		String aabits=null;
		try{
			header(in,"#bbgap_alphabet_config\t1");
			final byte[] line=in.nextLine();
			require(line!=null,"Compact alphabet resource must contain one aabits mapping");
			final String text=new String(line,StandardCharsets.US_ASCII);
			require(text.startsWith("aabits="),
				"Compact alphabet resource must contain one aabits mapping");
			aabits=text.substring(7);
			require(aabits.length()>0 && in.nextLine()==null,
				"Compact alphabet resource requires exactly one nonempty aabits line");
		}finally{require(!in.close(),"Alphabet reader failed");}
		return aabits;
	}
	private static void readScalarBounds(Path path,Expected e){
		final ByteFile in=ByteFile.makeByteFile(path.toString(),false);
		final LineParser1 lp=new LineParser1('\t');
		int rows=0;
		try{
			header(in,"scalar\trows\tmin\tmax\tmin_family\tmax_family");
			for(byte[] line=in.nextLine();line!=null;line=in.nextLine()){
				lp.set(line); require(lp.terms()==6,"Scalar bounds require six fields");
				final String name=lp.parseString(0);
				require(lp.parseLong(1)==79565,"Combined25 scalar bounds must have 79565 rows: "+name);
				require(rows<SCALAR_NAMES.length,"Too many combined25 scalar bound rows");
				require(name.equals(SCALAR_NAMES[rows]),"Unexpected scalar bound row: "+name);
				e.scalarBounds.put(name+"_min",lp.parseString(2));
				e.scalarBounds.put(name+"_max",lp.parseString(3));
				rows++;
			}
		}finally{require(!in.close(),"Scalar bounds reader failed");}
		require(rows==SCALAR_NAMES.length,"Scalar bounds must cover all combined25 rows");
	}
	private static char[][] readProfiles(String fname,FamilyReference[] refs){
		final HashMap<String,FamilyReference> byId=new HashMap<String,FamilyReference>(refs.length*2);
		for(FamilyReference f:refs){
			require(byId.put(f.id,f)==null,"Duplicate reference ID before profile read: "+f.id);
		}
		final char[][] profiles=new char[refs.length][];
		final ByteFile in=ByteFile.makeByteFile(fname,false);
		final LineParser1 lp=new LineParser1('\t');
		int count=0; long words=0;
		try{
			header(in,DDLFamilyProfiles.HEADER16);
			for(byte[] line=in.nextLine();line!=null;line=in.nextLine()){
				lp.set(line);
				require(lp.terms()==10 && lp.termEquals("tail",0) &&
					lp.parseInt(2)==FamilyCompactCodec.DDL_WIDTH,
					"Compact profile oracle requires 512-wide tail rows");
				require(lp.termEquals("literal",6) || lp.termEquals("nonempty",6),
					"Compact profile oracle row has unknown mode");
				if(lp.termEquals("literal",6)){continue;}
				final FamilyReference f=byId.get(lp.parseString(1));
				require(f!=null && f.members>0 && profiles[f.rank]==null &&
					lp.parseLong(3)==f.members && lp.parseLong(4)==f.residues &&
					lp.parseLong(5)==f.validWords,"Compact profile oracle identity/count mismatch");
				final char[] values=new char[FamilyCompactCodec.DDL_WIDTH];
				decimalCharArray(lp,7,values);
				profiles[f.rank]=values; count++; words=Math.addExact(words,f.validWords);
			}
		}finally{require(!in.close(),"Compact profile oracle read failed");}
		int expected=0; long expectedWords=0;
		for(FamilyReference f:refs){
			if(f.members>0){expected++; expectedWords=Math.addExact(expectedWords,f.validWords);}
		}
		require(count==expected && words==expectedWords,
			"Compact profile oracle must cover every live family exactly once");
		return profiles;
	}
	private static void header(ByteFile in,String expected){
		final byte[] line=in.nextLine();
		require(line!=null && expected.equals(new String(line,StandardCharsets.US_ASCII)),
			"Unexpected profile oracle header");
	}
	private static void decimalCharArray(LineParser1 lp,int field,char[] target){
		lp.setBounds(field);
		final byte[] bytes=lp.line(); int slot=0; long value=0; boolean digitSeen=false;
		for(int i=lp.a();i<=lp.b();i++){
			if(i==lp.b() || bytes[i]==','){
				require(digitSeen && slot<target.length && value<=65535,
					"Compact profile values must be 512 uint16 decimal terms");
				target[slot++]=(char)value; value=0; digitSeen=false;
			}else{
				final int digit=bytes[i]-'0';
				require(digit>=0 && digit<=9 && value<=(65535-digit)/10,
					"Invalid compact profile uint16 value");
				value=value*10+digit; digitSeen=true;
			}
		}
		require(slot==target.length,"Compact profile array must have 512 terms");
	}
	private static BlacklistFasta readBlacklist(Path path) throws Exception{
		final ArrayList<String> ids=new ArrayList<String>(), seqs=new ArrayList<String>();
		String id=null; final StringBuilder bases=new StringBuilder();
		final ByteFile in=ByteFile.makeByteFile(path.toString(),false);
		try{
			for(byte[] line=in.nextLine();line!=null;line=in.nextLine()){
				require(line.length>0,"Empty line in blacklist FASTA: "+path);
				if(line[0]=='>'){
					if(id!=null){ids.add(id); seqs.add(bases.toString());}
					id=new String(line,1,line.length-1,StandardCharsets.US_ASCII); bases.setLength(0);
				}else{
					require(id!=null,"Blacklist FASTA sequence before first header: "+path);
					bases.append(new String(line,StandardCharsets.US_ASCII));
				}
			}
		}finally{require(!in.close(),"Blacklist reader failed: "+path);}
		require(id!=null,"Blacklist FASTA has no records: "+path);
		ids.add(id); seqs.add(bases.toString());
		final BlacklistHeader header=new BlacklistHeader(ids.get(0));
		require(seqs.size()==header.words,"Header word count "+header.words+
			" does not match FASTA record count "+seqs.size()+": "+path);
		final long[] keys=new long[seqs.size()];
		for(int i=0;i<seqs.size();i++){
			final byte[] sequence=seqs.get(i).getBytes(StandardCharsets.US_ASCII);
			require(sequence.length==header.k,"Blacklist record "+(i+1)+" length "+sequence.length+
				" differs from header k="+header.k+": "+path);
			keys[i]=header.geometry.word(DDLPentamerBlacklist.pack(sequence,header.k));
		}
		Arrays.sort(keys);
		for(int i=1;i<keys.length;i++){
			require(keys[i-1]!=keys[i],"Duplicate reduced blacklist key "+keys[i]+": "+path);
		}
		DDLPentamerBlacklist.loadReducedStructural(path,"unmeasured",header.words,
			header.geometry,header.provenance);
		return new BlacklistFasta(path,ids,seqs,header,keys);
	}
	private static Counts writeNeutralBlacklist(BlacklistFasta a,BlacklistFasta b,Path out) throws Exception{
		require(a.header.samePolicy(b.header),"Blacklist headers differ beyond provenance: "+a.path+" vs "+b.path);
		require(a.seqs.size()==2048 && a.header.k==6 && "reduced_literal".equals(a.header.mode),
			"Q76 needs a K6 reduced_literal 2048-word blacklist");
		require(Arrays.equals(a.keys,b.keys),"SHORT and TAILS blacklist loaded keys differ");
		for(int i=0;i<a.seqs.size();i++){
			if(i>0){require(a.ids.get(i).equals(b.ids.get(i)),"Blacklist record ID "+(i+1)+" differs");}
			require(a.seqs.get(i).equals(b.seqs.get(i)),"Blacklist sequence "+(i+1)+" differs");
		}
		final String first="bbgap_blacklist_v1;k="+a.header.k+";words="+a.header.words+
			";mode="+a.header.mode+";provenance=bbgap_blacklist_k6_2048_20261116_v1"+
			";aabits="+a.header.aabitsRaw;
		try(BufferedWriter writer=Files.newBufferedWriter(out,StandardCharsets.US_ASCII)){
			for(int i=0;i<a.seqs.size();i++){
				writer.write('>'); writer.write(i==0 ? first : a.ids.get(i)); writer.newLine();
				writer.write(a.seqs.get(i)); writer.newLine();
			}
		}
		final BlacklistFasta neutral=readBlacklist(out);
		require("bbgap_blacklist_k6_2048_20261116_v1".equals(neutral.header.provenance),
			"Neutral blacklist provenance projection failed");
		require(Arrays.equals(a.keys,neutral.keys),"Neutral blacklist loaded keys differ");
		return new Counts(neutral.seqs.size(),neutral.seqs.size(),0,0,0);
	}
	private static void equal(String a,String b,String label,int ordinal){
		require(a.equals(b),label+" mismatch at ordinal "+ordinal+": "+a+" != "+b);
	}
	private static String scoreBinding(ProteinSubstitutionMatrix matrix){
		final String binding=matrix.binding();
		return binding.substring(binding.indexOf(':')+1);
	}
	private static int raw(float value){return Float.floatToRawIntBits(value);}
	private static long raw(double value){return Double.doubleToRawLongBits(value);}
	private static boolean rawEquals(float[] a,float[] b){
		if(a==null || b==null || a.length!=b.length){return false;}
		for(int i=0;i<a.length;i++){if(raw(a[i])!=raw(b[i])){return false;}}
		return true;
	}
	private static boolean rawEquals(double[] a,double[] b){
		if(a==null || b==null || a.length!=b.length){return false;}
		for(int i=0;i<a.length;i++){if(raw(a[i])!=raw(b[i])){return false;}}
		return true;
	}
	private static void require(boolean condition,String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	private static final class Sidecars{
		Sidecars(int size){
			lengthMedians=new double[size]; kdMeans=new double[size];
			gcMeans=new double[size]; gcSampleSds=new double[size]; stableIds=new int[size];
			seenLength=new boolean[size]; seenHydro=new boolean[size]; seenGc=new boolean[size];
		}
		final double[] lengthMedians,kdMeans,gcMeans,gcSampleSds;
		final int[] stableIds;
		final boolean[] seenLength,seenHydro,seenGc;
	}

	private static final class Counts{
		Counts(int physical_,int live_,long members_,long residues_,long validWords_){
			physical=physical_; live=live_; members=members_; residues=residues_; validWords=validWords_;
		}
		String tsv(){
			return "field\tvalue\nphysical\t"+physical+"\nlive\t"+live+"\nempty\t"+
				(physical-live)+"\nmembers\t"+members+"\nresidues\t"+residues+
				"\nvalid_words\t"+validWords+"\n";
		}
		final int physical,live;
		final long members,residues,validWords;
	}

	private static final class Expected{
		Expected(String view_,int families_,int live_,int empty_,long members_,int minAa_,
				int maxAa_,int alignmentMinSpan_,int alignmentMinPaired_,String marginScale_){
			view=view_; families=families_; live=live_; empty=empty_; members=members_;
			minAa=minAa_; maxAa=maxAa_; alignmentMinSpan=alignmentMinSpan_;
			alignmentMinPaired=alignmentMinPaired_; marginScale=marginScale_;
		}
		final String view,marginScale;
		final int families,live,empty,minAa,maxAa,alignmentMinSpan,alignmentMinPaired;
		final long members;
		ProteinKmerGeometry geometry;
		String matrixBinding;
		final HashMap<String,String> scalarBounds=new HashMap<String,String>();
	}

	private static final class BlacklistFasta{
		BlacklistFasta(Path path_,ArrayList<String> ids_,ArrayList<String> seqs_,
				BlacklistHeader header_,long[] keys_){
			path=path_; ids=ids_; seqs=seqs_; header=header_; keys=keys_;
		}
		final Path path;
		final ArrayList<String> ids,seqs;
		final BlacklistHeader header;
		final long[] keys;
	}

	private static final class BlacklistHeader{
		BlacklistHeader(String id){
			int k_=-1,words_=-1; String mode_=null,provenance_=null,aabits_=null,aabitsRaw_=null;
			int fields=0;
			for(int start=0;start<=id.length();){
				final int semi=id.indexOf(';',start);
				final int stop=semi<0 ? id.length() : semi;
				final String field=id.substring(start,stop);
				if(fields==0){
					require(field.equals("bbgap_blacklist_v1"),"Unknown BBGAP blacklist header: "+id);
				}else{
					final int equals=field.indexOf('=');
					require(equals>0 && equals<field.length()-1,
						"Malformed BBGAP blacklist header field: "+field);
					final String key=field.substring(0,equals), raw=field.substring(equals+1);
					final String value=decodeHeaderValue(raw);
					if("k".equals(key)){k_=Integer.parseInt(value);}
					else if("words".equals(key)){words_=Integer.parseInt(value);}
					else if("mode".equals(key)){mode_=value;}
					else if("provenance".equals(key)){provenance_=value;}
					else if("aabits".equals(key)){aabits_=value; aabitsRaw_=raw;}
					else{throw new IllegalArgumentException("Unknown BBGAP blacklist header field: "+key);}
				}
				fields++;
				if(semi<0){break;}
				start=semi+1;
			}
			require(fields>=5,"Incomplete BBGAP blacklist header: "+id);
			require(k_==6 && words_==2048 && "reduced_literal".equals(mode_) &&
				provenance_!=null && aabits_!=null,"Incomplete or unexpected blacklist header");
			k=k_; words=words_; mode=mode_; provenance=provenance_; aabits=aabits_; aabitsRaw=aabitsRaw_;
			geometry=ProteinKmerGeometry.make(k,aabits);
		}
		boolean samePolicy(BlacklistHeader b){
			return k==b.k && words==b.words && mode.equals(b.mode) && aabits.equals(b.aabits);
		}
		static String decodeHeaderValue(String value){
			final StringBuilder decoded=new StringBuilder(value.length());
			for(int i=0;i<value.length();i++){
				final char c=value.charAt(i);
				if(c=='%'){
					require(i+2<value.length(),"Truncated BBGAP blacklist header escape: "+value);
					final int hi=Character.digit(value.charAt(++i),16),lo=Character.digit(value.charAt(++i),16);
					require(hi>=0 && lo>=0,"Malformed BBGAP blacklist header escape: "+value);
					decoded.append((char)((hi<<4)|lo));
				}else{decoded.append(c);}
			}
			return decoded.toString();
		}
		final int k,words;
		final String mode,provenance,aabits,aabitsRaw;
		final ProteinKmerGeometry geometry;
	}

	private static final String[] SCALAR_NAMES={"hydrophobicity","net_charge","aromaticity","fraction_small","pI_tail"};
}
