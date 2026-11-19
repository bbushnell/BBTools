package test.ddl;

import java.io.ByteArrayOutputStream;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.nio.file.attribute.FileTime;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.Random;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import ddl.DDLPentamerBlacklist;
import ddl.ProteinKmerCandidates;
import ddl.ProteinKmerGeometry;
import dna.Data;
import fileIO.ReadWrite;

/** Literal selection oracles, independent FASTA rendering, malformed-input and stability fixtures.
 * Generated words are software fixtures, never biological calibration/evaluation populations.
 * @author Yaoyao
 */
public final class ProteinKmerCandidatesTest{
	public static void main(String[] args) throws Exception{
		check(args.length==1,"Require fresh fixture directory");final Path root=Paths.get(args[0]);
		check(root.isAbsolute() && !Files.exists(root),"Require new absolute fixture directory");Files.createDirectory(root);
		selection();ranked(root);stability(root);endToEnd(root);fullCounts(root);extended(root);
		ProteinKmerWeightedCandidatesTest.run(root.resolve("weighted"));weightedPackage(root);k5Package(root);
		System.out.println("PROTEIN_KMER_CANDIDATES_TEST_PASS checks="+checks);
	}

	private static void selection() throws Exception{
		equal(ProteinKmerCandidates.parseSizes("2049,1,40961"),new int[]{1,2048,2049,40961},"Arbitrary sizes sorted with retained baseline");
		equal(ProteinKmerCandidates.parseSizes("2048"),new int[]{2048},"Explicit baseline is not duplicated");
		for(String invalid:new String[]{"","0","-1","2049,2049","2049,",",2049","64000001","999999999999","1.5"}){
			reject(()->ProteinKmerCandidates.parseSizes(invalid),"Invalid size list "+invalid);
		}
		equal(ProteinKmerCandidates.mixed(new int[]{10,20,30},new int[]{10,40,50},4),new int[]{10,20,30,40},"Duplicate must consume its encounter");
		equal(ProteinKmerCandidates.mixed(new int[]{90,80},new int[]{10,20},1),new int[]{90},"Real first and immediate stop");
		equal(ProteinKmerCandidates.mixed(new int[]{90,80},new int[]{10,20},3),new int[]{10,80,90},"Emission order differs from selection order");
		equal(ProteinKmerCandidates.mixed(new int[]{1},new int[]{1,2,3},3),new int[]{1,2,3},"Exhausted real side must not truncate shadow");
		equal(ProteinKmerCandidates.mixed(new int[]{},new int[]{7,8},2),new int[]{7,8},"Empty side");
		reject(()->ProteinKmerCandidates.mixed(new int[]{1,2},new int[]{2,1},3),"Insufficient distinct candidates");
		reject(()->ProteinKmerCandidates.mixed(new int[]{},new int[]{},1),"Both exhausted");
		reject(()->ProteinKmerCandidates.mixed(new int[]{1},new int[]{2},0),"Zero limit");
		final int[] real=new int[2048],shadow=new int[2048],expected=new int[2048];
		for(int i=0;i<2048;i++){real[i]=shadow[i]=expected[i]=i;}
		shadow[2047]=9999;
		equal(ProteinKmerCandidates.mixed(real,shadow,2048),expected,"Stop on real encounter2048 before final shadow9999");
		for(int i=0;i<2048;i++){real[i]=2*i;shadow[i]=2*i+1;}
		equal(ProteinKmerCandidates.mixed(real,shadow,2048),expected,"Stop on shadow encounter at exact2048 boundary");
		final int[] duplicateReal=new int[40960],duplicateShadow=new int[40960];
		for(int i=0;i<40960;i++){duplicateReal[i]=duplicateShadow[i]=i;}duplicateShadow[40959]=999999;
		for(int n:new int[]{8192,40960}){equal(ProteinKmerCandidates.mixed(duplicateReal,duplicateShadow,n),Arrays.copyOf(duplicateReal,n),"Large nested duplicate-heavy mixed boundary");}
		for(int i=0;i<40960;i++){duplicateReal[i]=duplicateShadow[i]=40959-i;}duplicateShadow[40959]=999999;
		for(int n:new int[]{2048,8192,40960}){
			final int[] descendingExpected=new int[n];for(int i=0;i<n;i++){descendingExpected[i]=40960-n+i;}
			equal(ProteinKmerCandidates.mixed(duplicateReal,duplicateShadow,n),descendingExpected,"Descending duplicate-heavy encounters; sort only selected prefix, stop before final shadow");
		}
	}

	private static void ranked(Path root) throws Exception{
		final Path table=root.resolve("tiny.tsv"),fasta=root.resolve("tiny.faa");
		final String header="rank\tpacked_word\tword\toccurrences\n";
		write(table,header+"1\t1\tAAAAR\t9\n2\t0\tAAAAA\t8\n3\t2\tAAAAN\t8\n");
		write(fasta,">1_count_9\nAAAAR\n>2_count_8\nAAAAA\n>3_count_8\nAAAAN\n");
		equal(ProteinKmerCandidates.readRanked(table,fasta,5,3),new int[]{1,0,2},"Keep frequency order for selection");
		write(table,header+"1\t1\tAAAAR\t9\n2\t2\tAAAAN\t8\n3\t0\tAAAAA\t8\n");
		reject(()->ProteinKmerCandidates.readRanked(table,fasta,5,3),"Reversed tie order");
		write(table,header+"1\t1\tAAAAR\t9\n2\t1\tAAAAR\t8\n3\t0\tAAAAA\t7\n");
		reject(()->ProteinKmerCandidates.readRanked(table,fasta,5,3),"Duplicate with different counts");
		write(table,header+"1\t1\tAAAAR\t9\n2\t0\tAAAAA\t8\n3\t2\tAAAAN\t8\n");
		write(fasta,">1_count_9\nAAAAR\n>2_count_7\nAAAAA\n>3_count_8\nAAAAN\n");
		reject(()->ProteinKmerCandidates.readRanked(table,fasta,5,3),"Wrong FASTA count header");
		write(fasta,">1_count_9\nAAAAR\n>2_count_8\nAAAAA\n>3_count_8\nAAAAX\n");
		reject(()->ProteinKmerCandidates.readRanked(table,fasta,5,3),"Noncanonical FASTA");
		write(table,header+"1\t2\tAAAAR\t9\n");
		reject(()->ProteinKmerCandidates.readRanked(table,fasta,5,1),"Incorrect packed key");
		write(table,header+"1\t0\tAAAAA\t0\n");
		reject(()->ProteinKmerCandidates.readRanked(table,fasta,5,1),"Zero occurrences");
		write(table,header+"2\t0\tAAAAA\t1\n");
		reject(()->ProteinKmerCandidates.readRanked(table,fasta,5,1),"Wrong rank");
		write(table,header+"1\t658099827\tVVVVVV\t1\n");write(fasta,">1_count_1\nVVVVVV\n");
		equal(ProteinKmerCandidates.readRanked(table,fasta,6,1),new int[]{658099827},"Retain all30 K6 key bits");
	}

	private static void stability(Path root) throws Exception{
		final Path path=root.resolve("stable.txt");write(path,"abc");
		new ProteinKmerCandidates.StableInput(path).verify();checks++;
		final ProteinKmerCandidates.StableInput size=new ProteinKmerCandidates.StableInput(path);write(path,"abcd");
		reject(size::verify,"Changed size");
		final ProteinKmerCandidates.StableInput mtime=new ProteinKmerCandidates.StableInput(path);
		Files.setLastModifiedTime(path,FileTime.fromMillis(Files.getLastModifiedTime(path).toMillis()+10000));
		reject(mtime::verify,"Changed mtime");
		final ProteinKmerCandidates.StableInput key=new ProteinKmerCandidates.StableInput(path);
		final FileTime time=Files.getLastModifiedTime(path);Files.move(path,root.resolve("stable.old.txt"));write(path,"abcd");Files.setLastModifiedTime(path,time);
		reject(key::verify,"Changed inode with identical size/mtime");
		final Path link=root.resolve("input.link");Files.createSymbolicLink(link,path);
		reject(()->new ProteinKmerCandidates.StableInput(link),"Symlink input");
		reject(()->new ProteinKmerCandidates.StableInput(root),"Directory input");
	}

	private static void endToEnd(Path root) throws Exception{
		final Path real=root.resolve("real"),shadow=root.resolve("shadow"),legacy=root.resolve("legacy.faa"),out=root.resolve("candidates");
		Files.createDirectory(real);Files.createDirectory(shadow);
		for(int k:new int[]{5,6}){phase(real,"real",k,0);phase(shadow,"shadow",k,4096);}
		write(legacy,fasta(0,2048,5,1));
		final String[] args={"real="+real,"shadow="+shadow,"legacy="+legacy,"out="+out,"t=1"};
		ProteinKmerCandidates.main(args);
		check(text(out.resolve("PASS")).equals("PROTEIN_KMER_CANDIDATES_PASS lists=9 adopted=false\n"),"Missing success marker");
		for(int k:new int[]{5,6}){
			check(text(out.resolve("real.k"+k+".faa")).equals(fasta(0,2048,k,1)),"Pure real set");
			check(text(out.resolve("shadow.k"+k+".faa")).equals(fasta(4096,2048,k,1)),"Pure shadow set");
			check(text(out.resolve("mixed.k"+k+".faa")).equals(fasta(0,1024,k,1)+fasta(4096,1024,k,1025)),"Mixed exact independent set");
			check(text(out.resolve("union.k"+k+".faa")).equals(fasta(0,2048,k,1)+fasta(4096,2048,k,2049)),"Full union separate size");
		}
		check(text(out.resolve("legacy.k5.faa")).equals(fasta(0,2048,5,1)),"Legacy exact same word set");
		check(text(out.resolve("lists.tsv")).contains("mixed.k6.faa\t6\t2048\tmixed\t1024\t1024\t0\tunmeasured\n"),"Membership accounting");
		check(text(out.resolve("recipe.txt")).contains("file_hashes_computed=0\n") && text(out.resolve("recipe.txt")).contains("adopted=false\n"),"Identity/adoption declarations");
		reject(()->ProteinKmerCandidates.main(args),"Output overwrite");
		final Path phase=real.resolve("real_k5"),result=phase.resolve("result");
		final String sizes=text(phase.resolve("accepted.output_sizes.tsv"));
		for(String name:new String[]{"recipe.txt","top.tsv.gz","top.faa.gz"}){
			final Path consumed=result.resolve(name);final byte[] original=Files.readAllBytes(consumed);
			if(name.equals("recipe.txt")){write(consumed,text(consumed)+"fixture_note=valid_extra_recipe_field\n");}
			else{try(GZIPOutputStream extra=new GZIPOutputStream(Files.newOutputStream(consumed,StandardOpenOption.APPEND))){extra.finish();}}
			check(Files.size(consumed)!=original.length,"Replacement must change stored size");
			check(text(phase.resolve("accepted.output_sizes.tsv")).equals(sizes) && text(phase.resolve("output_sizes.tsv")).equals(sizes),"Receipts must stay unchanged");
			final Path rejected=root.resolve("replaced_"+name);final String[] replacement=args.clone();replacement[3]="out="+rejected;
			rejectMessage(()->ProteinKmerCandidates.main(replacement),"Consumed file differs from accepted size:");
			check(!Files.exists(rejected),"Changed accepted file published output");Files.write(consumed,original);
		}
		final String[] invalidReceipts={sizes.replace("file\tbytes", "file\tsize"),
			sizes.replace("top.faa.gz\t", "top.tsv.gz\t"),
			sizes.replace("recipe.txt\t"+Files.size(result.resolve("recipe.txt"))+"\n", "")};
		for(int i=0;i<invalidReceipts.length;i++){
			write(phase.resolve("accepted.output_sizes.tsv"),invalidReceipts[i]);write(phase.resolve("output_sizes.tsv"),invalidReceipts[i]);
			final Path rejected=root.resolve("invalid_size_receipt_"+i);final String[] invalid=args.clone();invalid[3]="out="+rejected;
			reject(()->ProteinKmerCandidates.main(invalid),"Invalid size receipt must fail even when saved copies agree");
			check(!Files.exists(rejected),"Invalid size receipt published output");
		}
		write(phase.resolve("accepted.output_sizes.tsv"),sizes);write(phase.resolve("output_sizes.tsv"),sizes);
		write(real.resolve("real_k5/accepted.time.tsv"),"changed\n");
		final String[] bad=args.clone();bad[3]="out="+root.resolve("bad_receipt");
		reject(()->ProteinKmerCandidates.main(bad),"Changed owner receipt");
		check(!Files.exists(root.resolve("bad_receipt")),"Rejected receipt published output");
	}


	/** Small independent full-sort oracle, then a literal large-boundary oracle. */
	private static void fullCounts(Path root) throws Exception{
		final Path table=root.resolve("count_oracle.tsv");final Random random=new Random(613);
		for(int k:new int[]{5,6}){
			final long[] counts=new long[257];final ArrayList<Integer> expected=new ArrayList<Integer>();
			final StringBuilder rows=new StringBuilder("packed_word\tword\toccurrences\n");long total=0;
			for(int i=0;i<counts.length;i++){
				counts[i]=1+random.nextInt(11);total+=counts[i];expected.add(i);final String word=word(i,k);
				rows.append(pack(word)).append('\t').append(word).append('\t').append(counts[i]).append('\n');
			}
			write(table,rows.toString());expected.sort((a,b)->counts[a]==counts[b] ? Integer.compare(a,b) : Long.compare(counts[b],counts[a]));
			for(int limit:new int[]{1,2,17,256,257}){
				final int[] actual=ProteinKmerCandidates.rankCounts(table,k,limit,257,total),oracle=new int[limit];
				for(int i=0;i<limit;i++){oracle[i]=pack(word(expected.get(i),k));}equal(actual,oracle,"Bounded heap/full-sort oracle K"+k+" limit"+limit);
			}
		}
		final String header="packed_word\tword\toccurrences\n";
		write(table,header+"0\tAAAAAA\t3\n1\tAAAAAR\t2\n2\tAAAAAN\t2\n");
		equal(ProteinKmerCandidates.rankCounts(table,6,2,3,7),new int[]{0,1},"Cutoff tie retains lower packed key");
		reject(()->ProteinKmerCandidates.rankCounts(table,6,2,4,7),"Missing count row");
		reject(()->ProteinKmerCandidates.rankCounts(table,6,2,3,8),"Incorrect accepted occurrence sum");
		reject(()->ProteinKmerCandidates.rankCounts(table,6,2,2,7),"Extra count row");
		reject(()->ProteinKmerCandidates.rankCounts(table,6,64000001,64000001,64000001),"K6 universe bound");
		write(table,header+"0\tAAAAAA\t2\n0\tAAAAAA\t2\n");
		reject(()->ProteinKmerCandidates.rankCounts(table,6,1,2,4),"Duplicate full-table key");
		write(table,header+"1\tAAAAAR\t2\n0\tAAAAAA\t2\n");
		reject(()->ProteinKmerCandidates.rankCounts(table,6,1,2,4),"Decreasing full-table key");
		write(table,header+"1\tAAAAAA\t2\n");
		reject(()->ProteinKmerCandidates.rankCounts(table,6,1,1,2),"Mismatched full-table key");
		write(table,header+"0\tAAAAAX\t2\n");
		reject(()->ProteinKmerCandidates.rankCounts(table,6,1,1,2),"Noncanonical full-table word");
		write(table,header+"0\tAAAAAA\t0\n");
		reject(()->ProteinKmerCandidates.rankCounts(table,6,1,1,1),"Nonpositive full-table count");
		write(table,header+"658099827\tVVVVVV\t1\n");
		equal(ProteinKmerCandidates.rankCounts(table,6,1,1,1),new int[]{658099827},"Full-table K6 high key");
		write(table,header+"0\tAAAAAA\t9223372036854775806\n1\tAAAAAR\t2\n");
		try{ProteinKmerCandidates.rankCounts(table,6,1,2,Long.MAX_VALUE);throw new AssertionError("Occurrence overflow must fail");}
		catch(ArithmeticException expected){checks++;}
		final StringBuilder large=new StringBuilder(header);final int[] boundary=new int[40960];
		for(int i=0;i<41000;i++){final String word=word(i,6);large.append(pack(word)).append('\t').append(word).append('\t').append(i+1).append('\n');}
		write(table,large.toString());
		for(int i=0;i<boundary.length;i++){boundary[i]=pack(word(40999-i,6));}
		equal(ProteinKmerCandidates.rankCounts(table,6,40960,41000,41000L*41001/2),boundary,"40960 heap boundary and repeated root replacement");
		final int[] above=Arrays.copyOf(boundary,40961);above[40960]=pack(word(39,6));
		equal(ProteinKmerCandidates.rankCounts(table,6,40961,41000,41000L*41001/2),above,"Selection above former hardcoded capacity");
	}

	/** Full counts must reproduce the nine default lists and all larger literal sets. */
	private static void extended(Path root) throws Exception{
		final Path inputs=root.resolve("extended_inputs"),real=inputs.resolve("real"),shadow=inputs.resolve("shadow");
		Files.createDirectory(inputs);Files.createDirectory(real);Files.createDirectory(shadow);
		phase(real,"real",5,0);phase(shadow,"shadow",5,4096);
		phase(real,"real",6,0,40960);phase(shadow,"shadow",6,65536,40960);
		final Path legacy=inputs.resolve("legacy.faa"),small=root.resolve("default_from_full"),large=root.resolve("extended_output");
		write(legacy,fasta(0,2048,5,1));
		final String[] original={"real="+real,"shadow="+shadow,"legacy="+legacy,"out="+small,"t=1"};
		ProteinKmerCandidates.main(original);
		final String[] args=Arrays.copyOf(original,6);args[3]="out="+large;args[5]="k6sizes=2048,8192,40960";
		check(Data.BGZIP() && Data.BGZIP_VERSION_levelFlag && Data.BGZIP_VERSION_threadsFlag,"Fixture requires bgzip with explicit level/thread flags");
		final boolean levelFlag=Data.BGZIP_VERSION_levelFlag,threadsFlag=Data.BGZIP_VERSION_threadsFlag;
		final String[] unavailable=args.clone();unavailable[3]="out="+root.resolve("unsupported_bgzip");
		try{
			Data.BGZIP_VERSION_levelFlag=false;
			rejectMessage(()->ProteinKmerCandidates.main(unavailable),"Extended candidates require external bgzip");
			Data.BGZIP_VERSION_levelFlag=levelFlag;Data.BGZIP_VERSION_threadsFlag=false;
			rejectMessage(()->ProteinKmerCandidates.main(unavailable),"Extended candidates require external bgzip");
		}finally{Data.BGZIP_VERSION_levelFlag=levelFlag;Data.BGZIP_VERSION_threadsFlag=threadsFlag;}
		check(!Files.exists(root.resolve("unsupported_bgzip")),"Missing backend capabilities must fail before publishing output");
		// Conflicting caller preferences must not weaken the selected backend or leak changed settings.
		final boolean bgzf=ReadWrite.USE_BGZF,bgzip=ReadWrite.USE_BGZIP,force=ReadWrite.FORCE_BGZIP,
			nativeBgzf=ReadWrite.ALLOW_NATIVE_BGZF,preferNative=ReadWrite.PREFER_NATIVE_BGZF_OUT,change=ReadWrite.ALLOW_ZIPLEVEL_CHANGE;
		final int level=ReadWrite.ZIPLEVEL;
		try{
			ReadWrite.USE_BGZF=ReadWrite.USE_BGZIP=ReadWrite.FORCE_BGZIP=false;
			ReadWrite.ALLOW_NATIVE_BGZF=ReadWrite.PREFER_NATIVE_BGZF_OUT=ReadWrite.ALLOW_ZIPLEVEL_CHANGE=true;ReadWrite.ZIPLEVEL=1;
			ProteinKmerCandidates.main(args);
			check(!ReadWrite.USE_BGZF && !ReadWrite.USE_BGZIP && !ReadWrite.FORCE_BGZIP && ReadWrite.ALLOW_NATIVE_BGZF &&
				ReadWrite.PREFER_NATIVE_BGZF_OUT && ReadWrite.ALLOW_ZIPLEVEL_CHANGE && ReadWrite.ZIPLEVEL==1,"Restore caller compression preferences");
		}finally{
			ReadWrite.USE_BGZF=bgzf;ReadWrite.USE_BGZIP=bgzip;ReadWrite.FORCE_BGZIP=force;
			ReadWrite.ALLOW_NATIVE_BGZF=nativeBgzf;ReadWrite.PREFER_NATIVE_BGZF_OUT=preferNative;
			ReadWrite.ALLOW_ZIPLEVEL_CHANGE=change;ReadWrite.ZIPLEVEL=level;
		}
		check(text(large.resolve("PASS")).equals("PROTEIN_KMER_CANDIDATES_PASS lists=15 adopted=false\n"),"Missing expanded success marker");
		for(String name:new String[]{"real.k5","shadow.k5","mixed.k5","union.k5","real.k6","shadow.k6","mixed.k6","union.k6","legacy.k5"}){
			check(Arrays.equals(Files.readAllBytes(small.resolve(name+".faa")),Files.readAllBytes(large.resolve(name+".faa"))),"Exact default2048/union/legacy parity: "+name);
		}
		for(int n:new int[]{8192,40960}){
			final String suffix=n==40960 ? ".k6.faa.gz" : ".k6.faa";
			check(decoded(large.resolve("real_"+n+suffix)).equals(fasta(0,n,6,1)),"Independent larger pure real set");
			check(decoded(large.resolve("shadow_"+n+suffix)).equals(fasta(65536,n,6,1)),"Independent larger pure shadow set");
			check(decoded(large.resolve("mixed_"+n+suffix)).equals(fasta(0,n/2,6,1)+fasta(65536,n/2,6,n/2+1)),"Independent larger mixed set and emission order");
		}
		for(String population:new String[]{"real","shadow","mixed"}){
			final byte[] compressed=Files.readAllBytes(large.resolve(population+"_40960.k6.faa.gz"));
			check(compressed.length>=18 && compressed[0]==31 && (compressed[1]&255)==139 && (compressed[3]&4)!=0 && compressed[12]=='B' && compressed[13]=='C',"Large lists require BGZF");
			final String expected=population.equals("real") ? fasta(0,40960,6,1) : population.equals("shadow") ? fasta(65536,40960,6,1) :
				fasta(0,20480,6,1)+fasta(65536,20480,6,20481);
			final Path oracle=root.resolve(population+"_level9_oracle.gz"),stderr=root.resolve(population+"_level9_oracle.stderr");
			// Feed independent literal bytes to stdin, matching the production pipe without using its writer/dispatch.
			final Process compressor=new ProcessBuilder("bgzip","-c","-f","-l","9","-@","1")
				.redirectOutput(oracle.toFile()).redirectError(stderr.toFile()).start();
			try(java.io.OutputStream stdin=compressor.getOutputStream()){stdin.write(expected.getBytes(StandardCharsets.US_ASCII));}
			check(compressor.waitFor()==0,"Independent explicit level9 compressor failed: "+stderr);
			check(Arrays.equals(compressed,Files.readAllBytes(oracle)),"Exact explicit bgzip -l9 byte oracle: "+population);
		}
		check(text(large.resolve("recipe.txt")).contains("full_real_top2048_parity=exact\n") &&
			text(large.resolve("recipe.txt")).contains("full_shadow_top2048_parity=exact\n"),"Both count populations require exact top parity");
		final Path odd=root.resolve("arbitrary_sizes");final String[] arbitrary=args.clone();
		arbitrary[3]="out="+odd;arbitrary[5]="k6sizes=2049,1";ProteinKmerCandidates.main(arbitrary);
		check(text(odd.resolve("PASS")).equals("PROTEIN_KMER_CANDIDATES_PASS lists=15 adopted=false\n"),"Arbitrary-size run completes");
		check(text(odd.resolve("recipe.txt")).contains("\nk6sizes=1,2048,2049\n") &&
			text(odd.resolve("recipe.txt")).contains("\nfull_count_selection=bounded_primitive_heap_top2049\n"),"Canonical sizes and exact count-read capacity");
		for(String name:new String[]{"real.k5","shadow.k5","mixed.k5","union.k5","real.k6","shadow.k6","mixed.k6","union.k6","legacy.k5"}){
			check(Arrays.equals(Files.readAllBytes(small.resolve(name+".faa")),Files.readAllBytes(odd.resolve(name+".faa"))),"Arbitrary size preserves default bytes: "+name);
		}
		check(text(odd.resolve("real_2049.k6.faa")).equals(fasta(0,2049,6,1)),"Exact2049 real words");
		check(text(odd.resolve("shadow_2049.k6.faa")).equals(fasta(65536,2049,6,1)),"Exact2049 shadow words");
		check(text(odd.resolve("mixed_2049.k6.faa")).equals(fasta(0,1025,6,1)+fasta(65536,1024,6,1026)),"Odd mixed size takes1025 real and1024 shadow");
		check(text(odd.resolve("mixed_1.k6.faa")).equals(fasta(0,1,6,1)),"Below-baseline size stops on first real word");
		final Path phase=real.resolve("real_k6"),table=phase.resolve("result/counts.tsv.gz");
		final byte[] saved=Files.readAllBytes(table);
		try(GZIPOutputStream extra=new GZIPOutputStream(Files.newOutputStream(table,StandardOpenOption.APPEND))){extra.finish();}
		final String[] bad=args.clone();bad[3]="out="+root.resolve("changed_full_size");
		rejectMessage(()->ProteinKmerCandidates.main(bad),"Consumed file differs from accepted size:");Files.write(table,saved);
		// Same population totals and syntactically valid ranking, but different accepted top counts.
		final Path top=phase.resolve("result/top.tsv.gz"),topFasta=phase.resolve("result/top.faa.gz");
		gzip(top,decoded(top).replace("\t2\n","\t3\n"));gzip(topFasta,decoded(topFasta).replace("_count_2\n","_count_3\n"));receipts(phase);
		final String[] mismatch=args.clone();mismatch[3]="out="+root.resolve("top_parity_mismatch");
		rejectMessage(()->ProteinKmerCandidates.main(mismatch),"Full-count/top2048 parity failure");
		check(!Files.exists(root.resolve("top_parity_mismatch")),"Failed parity must not publish candidates");
	}
	private static void weightedPackage(Path root)throws Exception{
		final Path real=root.resolve("weighted_real"),shadow=root.resolve("weighted_shadow"),out=root.resolve("weighted_package");
		Files.createDirectory(real);Files.createDirectory(shadow);phase(real,"real",6,0,4096);phase(shadow,"shadow",6,4096,2048);
		final String[] args={"real="+real,"shadow="+shadow,"out="+out,"arms=raw:shadow:2049,raw:weighted:2049,reduced:shadow:17,reduced:weighted:19",
			"alphabet=all:A/C/DE/FWY/G/IV/LM/NH/P/Q/RK/ST","negativeweight=4","positiveweight=1","t=1"};
		// Only2048 negative raw words exist: request must fail rather than pad or emit a shorter list.
		reject(()->ProteinKmerCandidates.main(args),"Insufficient raw negatives");check(!Files.exists(out),"Failed selector cannot publish package");
		args[3]="arms=raw:shadow:2048,raw:weighted:2049,reduced:shadow:17,reduced:weighted:19";
		// Separate treatment capacities permit weighted union sizes larger than negative-only support.
		ProteinKmerCandidates.main(args);
		check(text(out.resolve("PASS")).equals("PROTEIN_KMER_WEIGHTED_CANDIDATES_PASS lists=4 adopted=false\n"),"Weighted package terminal marker");
		check(text(out.resolve("raw_shadow_k6_2048.faa")).equals(fasta(4096,2048,6,1)),"Raw negative-only list preserves accepted words");
		check(text(out.resolve("raw_weighted_k6_2049.faa")).equals(fasta(0,1,6,1)+fasta(4096,2048,6,2)),"Normalized weighting selects negatives then lowest positive tie");
		final Path scorePath=out.resolve("raw_weighted_k6_2049.scores.tsv.gz");
		final String[] scoreRows=decoded(scorePath).split("\n",-1);
		check(scoreRows.length==2051 && scoreRows[0].equals("rank\tpacked_word\tnegative_occurrences\tpositive_occurrences\tinteger_score") &&
			scoreRows[1].equals("1\t"+pack(word(4096,6))+"\t2\t0\t16") && scoreRows[2049].equals("2049\t0\t0\t2\t2") && scoreRows[2050].isEmpty(),
			"Decoded score provenance preserves exact normalized ordering, both counts and integer coefficients");
		final byte[] scoreBytes=Files.readAllBytes(scorePath);
		check(scoreBytes.length>=18 && scoreBytes[0]==31 && (scoreBytes[1]&255)==139 && (scoreBytes[3]&4)!=0 && scoreBytes[12]=='B' && scoreBytes[13]=='C',
			"Score tables use BGZF directly");
		check(text(out.resolve("lists.tsv")).contains("\tweighted\treduced\tall:A/C/DE/FWY/G/IV/LM/NH/P/Q/RK/ST\t4\t1\t4096\t8192\tunmeasured\n"),"Reduced manifest records weights, alphabet and both denominators");
		check(text(out.resolve("recipe.txt")).contains("full_count_top2048_parity=exact_both_populations"),"Full count parity retained");
		final String[] versioned=Arrays.copyOf(args,args.length+2);versioned[2]="out="+root.resolve("weighted_versioned");
		versioned[versioned.length-2]="versionedblacklistheaders=t";versioned[versioned.length-1]="blacklistprovenance=protein_kmer_candidates";
		ProteinKmerCandidates.main(versioned);
		final Path reduced=Paths.get(versioned[2].substring(4)).resolve("reduced_shadow_k6_17.faa");
		check(text(reduced).startsWith(">bbgap_blacklist_v1;k=6;words=17;mode=reduced_literal;provenance=protein_kmer_candidates;aabits=all:A/C/DE/FWY/G/IV/LM/NH/P/Q/RK/ST\n"),
			"Versioned blacklist header carries K, count, mode, provenance and exact alphabet");
		final ProteinKmerGeometry geometry=ProteinKmerGeometry.make(6,"all:A/C/DE/FWY/G/IV/LM/NH/P/Q/RK/ST");
		check(DDLPentamerBlacklist.loadReducedStructural(reduced,"unmeasured",17,geometry,"protein_kmer_candidates").size()==17,
			"Versioned reduced blacklist remains loadable with the exact geometry");
		rankedExport(root,out);
		reject(()->ProteinKmerCandidates.main(args),"Output overwrite");
	}
	private static void k5Package(Path root)throws Exception{
		final Path real=root.resolve("k5_real"),shadow=root.resolve("k5_shadow"),out=root.resolve("k5_package");
		Files.createDirectory(real);Files.createDirectory(shadow);phase(real,"real",5,0,4096);phase(shadow,"shadow",5,8192,4096);
		final String[] args={"real="+real,"shadow="+shadow,"out="+out,"k=5",
			"arms=raw:real:1024,raw:shadow:4096,reduced:real:17,reduced:shadow:19",
			"alphabet=all:A/C/D/E/FWY/G/H/IV/K/LM/N/P/Q/R/S/T","t=1"};
		ProteinKmerCandidates.main(args);
		check(text(out.resolve("raw_real_k5_1024.faa")).equals(fasta(0,1024,5,1)),"K5 real counts selected with packed-key tie break");
		check(text(out.resolve("raw_shadow_k5_4096.faa")).equals(fasta(8192,4096,5,1)),"K5 shadow full counts selected beyond retained top2048");
		check(text(out.resolve("lists.tsv")).contains("raw_real_k5_1024.faa\t5\t1024\treal\traw\tnone\t0\t1\t8192\t8192\tunmeasured\n"),"Real-only manifest cannot claim negative weight");
		check(decoded(out.resolve("raw_real_k5_1024.scores.tsv.gz")).startsWith("rank\tpacked_word\tnegative_occurrences\tpositive_occurrences\tinteger_score\n1\t0\t0\t2\t2\n"),"Real-only score provenance");
		check(text(out.resolve("recipe.txt")).contains("\nk=5\n") && text(out.resolve("PASS")).contains("lists=4"),"K5 binds recipe and exact output count");
		final String[] bad=args.clone();bad[2]="out="+root.resolve("k5_bad_k");bad[3]="k=4";
		reject(()->ProteinKmerCandidates.main(bad),"Unsupported explicit K");check(!Files.exists(root.resolve("k5_bad_k")),"Invalid K publishes no output");
		final String[] defaultK=args.clone();defaultK[2]="out="+root.resolve("k5_as_default_k6");defaultK[3]="negativeweight=4";
		boolean missing=false;try{ProteinKmerCandidates.main(defaultK);}catch(java.nio.file.NoSuchFileException expected){missing=true;}
		check(missing && !Files.exists(root.resolve("k5_as_default_k6")),"Omitted K still requests K6 count inputs");
		reject(()->new ProteinKmerCandidates(new String[]{"rankedexport="+out,"out="+root.resolve("bad_export_k"),"k=5"}),"Ranked export cannot change K");
		reject(()->new ProteinKmerCandidates(new String[]{"real="+real,"shadow="+shadow,"legacy="+out,"out="+root.resolve("bad_old_k"),"k=5"}),"Legacy selector cannot silently switch K");
		reject(()->new ProteinKmerCandidates(new String[]{"rankedexport="+out,"out="+root.resolve("bad_export_header"),
			"versionedblacklistheaders=t"}),"Ranked export cannot accept unused blacklist headers");
		reject(()->new ProteinKmerCandidates(new String[]{"real="+real,"shadow="+shadow,"legacy="+out,
			"out="+root.resolve("bad_old_header"),"versionedblacklistheaders=t"}),"Legacy selector cannot accept unused blacklist headers");
	}

	/** Check rank-preserving export without count files, and reject inconsistent accepted artifacts. */
	private static void rankedExport(Path root,Path generated)throws Exception{
		final Path accepted=root.resolve("accepted_export_fixture"),result=accepted.resolve("result"),out=root.resolve("ranked_export");
		Files.createDirectory(accepted);Files.createDirectory(result);
		for(String name:new String[]{"lists.tsv","PASS","raw_shadow_k6_2048.faa","raw_weighted_k6_2049.faa",
			"raw_shadow_k6_2048.scores.tsv.gz","raw_weighted_k6_2049.scores.tsv.gz"}){Files.copy(generated.resolve(name),result.resolve(name));}
		write(accepted.resolve("accepted.utc"),"synthetic software fixture\n");
		write(accepted.resolve("native.exit"),"0\n");write(accepted.resolve("exit"),"0\n");
		write(accepted.resolve("PASS"),"PROTEIN_KMER_CANDIDATES_RUN_PASS lists=4 adopted=false\n");
		final String[] args={"rankedexport="+accepted,"out="+out,"t=1"};
		ProteinKmerCandidates.main(args);
		for(String selection:new String[]{"shadow","weighted"}){
			final int size=selection.equals("shadow") ? 2048 : 2049;final StringBuilder expected=new StringBuilder();
			for(int i=0;i<size;i++){
				final String word=word(i<2048 ? 4096+i : 0,6);
				expected.append(">rank_").append(i+1).append(" packed_word=").append(pack(word)).append('\n').append(word).append('\n');
			}
			check(decoded(out.resolve("raw_"+selection+"_k6_"+size+".ranked.faa.gz")).equals(expected.toString()),
				"Export preserves independently expected ranking, rather than original packed-word FASTA order: "+selection);
		}
		check(text(out.resolve("PASS")).equals("PROTEIN_KMER_RANKED_EXPORT_PASS lists=2 adopted=false\n"),"Export completes without source count files");
		reject(()->ProteinKmerCandidates.main(args),"Ranked export cannot overwrite output");
		final String[] ambiguous={args[0],"out="+root.resolve("ambiguous_export"),"arms=raw:shadow:1"};
		reject(()->ProteinKmerCandidates.main(ambiguous),"Export cannot silently reselect");
		final Path score=result.resolve("raw_shadow_k6_2048.scores.tsv.gz");final String valid=decoded(score);
		final String[] invalid={valid.replaceFirst("1\\t", "2\t"),
			valid.replaceFirst("\\t2\\t0\\t2\\n", "\t2\t0\t3\n"),valid.substring(0,valid.lastIndexOf('\n',valid.length()-2)+1)};
		for(int i=0;i<invalid.length;i++){
			gzip(score,invalid[i]);final Path bad=root.resolve("bad_export_score_"+i);
			reject(()->ProteinKmerCandidates.main(new String[]{args[0],"out="+bad}),"Reject bad rank, score or truncated ranking");
			check(!Files.exists(bad),"Invalid score table cannot publish ranked output");
		}gzip(score,valid);
		final Path fasta=result.resolve("raw_shadow_k6_2048.faa");final String original=text(fasta);
		write(fasta,original.replaceFirst(word(4096,6),"VVVVVV"));
		reject(()->ProteinKmerCandidates.main(new String[]{args[0],"out="+root.resolve("bad_export_membership")}),"Score and FASTA word sets must agree");
		write(fasta,original);
		final String manifest=text(result.resolve("lists.tsv"));
		write(result.resolve("lists.tsv"),manifest.replace("raw_shadow_k6_2048.faa","../raw_shadow_k6_2048.faa"));
		reject(()->ProteinKmerCandidates.main(new String[]{args[0],"out="+root.resolve("bad_export_path")}),"Reject manifest path traversal");
		write(result.resolve("lists.tsv"),manifest);
	}
	private static int pack(String word){int packed=0;for(int i=0;i<word.length();i++){packed=packed*32+AA.indexOf(word.charAt(i));}return packed;}
	private static String decoded(Path path) throws Exception{
		if(!path.toString().endsWith(".gz")){return text(path);}
		try(InputStream in=new GZIPInputStream(Files.newInputStream(path));ByteArrayOutputStream out=new ByteArrayOutputStream()){
			final byte[] buffer=new byte[8192];for(int n;(n=in.read(buffer))>=0;){out.write(buffer,0,n);}return new String(out.toByteArray(),StandardCharsets.US_ASCII);
		}
	}

	private static void phase(Path root,String population,int k,int first) throws Exception{phase(root,population,k,first,2048);}
	private static void phase(Path root,String population,int k,int first,int distinct) throws Exception{
		final Path p=root.resolve(population+"_k"+k),r=p.resolve("result");Files.createDirectory(p);Files.createDirectory(r);
		write(p.resolve("accepted.utc"),"software fixture\n");
		for(String f:new String[]{"native.exit","exit","watch.exit"}){write(p.resolve(f),"0\n");}
		write(p.resolve("PASS"),"PROTEIN_KMER_COUNTS_PASS\n");write(r.resolve("PASS"),"PROTEIN_KMER_FREQUENCY_PASS\n");
		for(String f:new String[]{"native.time.tsv","accepted.time.tsv"}){write(p.resolve(f),"fixture\n");}
		write(r.resolve("recipe.txt"),"in=software_fixture\npopulation="+population+"\nk="+k+"\ntop=2048\n");
		final StringBuilder table=new StringBuilder("rank\tpacked_word\tword\toccurrences\n"),fasta=new StringBuilder();
		for(int i=0;i<2048;i++){
			final String word=word(first+i,k);int packed=0;
			for(int j=0;j<word.length();j++){packed=packed*32+AA.indexOf(word.charAt(j));}
			table.append(i+1).append('\t').append(packed).append('\t').append(word).append("\t2\n");
			fasta.append('>').append(i+1).append("_count_2\n").append(word).append('\n');
		}
		gzip(r.resolve("top.tsv.gz"),table.toString());gzip(r.resolve("top.faa.gz"),fasta.toString());
		final StringBuilder counts=new StringBuilder("packed_word\tword\toccurrences\n");
		for(int i=0;i<distinct;i++){final String word=word(first+i,k);counts.append(pack(word)).append('\t').append(word).append("\t2\n");}
		gzip(r.resolve("counts.tsv.gz"),counts.toString());
		write(r.resolve("summary.tsv"),"population\tk\trecords\tresidues\tvalid_words\tdistinct_words\tselected_words\n"+
			population+"\t"+k+"\t1\t"+(2L*distinct+k-1)+"\t"+(2L*distinct)+"\t"+distinct+"\t2048\n");
		receipts(p);
	}
	private static void receipts(Path p) throws Exception{
		final Path r=p.resolve("result");
		final StringBuilder sizes=new StringBuilder("file\tbytes\n");
		for(String name:new String[]{"counts.tsv.gz","top.tsv.gz","top.faa.gz","summary.tsv","recipe.txt","PASS"}){
			sizes.append(name).append('\t').append(Files.size(r.resolve(name))).append('\n');
		}
		write(p.resolve("output_sizes.tsv"),sizes.toString());write(p.resolve("accepted.output_sizes.tsv"),sizes.toString());
	}
	private static String word(int index,int k){
		final char[] letters=new char[k];for(int i=k-1;i>=0;i--){letters[i]=AA.charAt(index%20);index/=20;}return new String(letters);
	}
	private static String fasta(int first,int count,int k,int firstHeader){
		final StringBuilder out=new StringBuilder();for(int i=0;i<count;i++){out.append('>').append(firstHeader+i).append('\n').append(word(first+i,k)).append('\n');}return out.toString();
	}
	private static void gzip(Path p,String s) throws Exception{try(GZIPOutputStream out=new GZIPOutputStream(Files.newOutputStream(p))){out.write(s.getBytes(StandardCharsets.US_ASCII));}}
	private static void write(Path p,String s) throws Exception{Files.write(p,s.getBytes(StandardCharsets.US_ASCII));}
	private static String text(Path p) throws Exception{return new String(Files.readAllBytes(p),StandardCharsets.US_ASCII);}
	private static void equal(int[] actual,int[] expected,String why){check(Arrays.equals(actual,expected),why);}
	private static void check(boolean ok,String why){if(!ok){throw new AssertionError(why);}checks++;}
	private static void reject(Action action,String why) throws Exception{
		try{action.run();}catch(IllegalArgumentException expected){checks++;return;}throw new AssertionError(why);
	}
	private static void rejectMessage(Action action,String prefix) throws Exception{
		try{action.run();}catch(IllegalArgumentException expected){check(expected.getMessage().startsWith(prefix),"Wrong rejection: "+expected);return;}
		throw new AssertionError("Expected rejection: "+prefix);
	}
	private interface Action{void run() throws Exception;}
	private static int checks;
	private static final String AA="ARNDCEQGHILKMFPSTWYV";
}
