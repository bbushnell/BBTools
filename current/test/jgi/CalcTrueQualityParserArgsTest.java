package test.jgi;

import java.io.IOException;
import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.ReadWrite;
import jgi.CalcTrueQuality;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import stream.SamLine;
import tracker.ReadStats;

/**
 * Constructor-only fixture for CalcTrueQuality argument parsing.
 * It creates readable sentinel inputs, instantiates CalcTrueQuality once per
 * case, and reflects parsed constructor fields without calling process(),
 * loading variants, starting workers, or writing calibration outputs.
 * @author Jean
 */
public final class CalcTrueQualityParserArgsTest{

	/** Runs all parser cases in the requested baseline or fixed phase. */
	public static void main(String[] args){
		final PreParser pp=new PreParser(args, CalcTrueQualityParserArgsTest.class, false);
		final SavedState initial=new SavedState();
		TestFiles files=null;
		try{
			final Args parsed=new Args(pp.args);
			files=TestFiles.create();
			new CalcTrueQualityParserArgsTest(parsed.fixed, files, initial).process();
		}finally{
			try{
				if(files!=null){files.cleanup();}
			}finally{
				try{
					initial.restoreAll();
				}finally{
					Shared.closeStream(pp.outstream);
				}
			}
		}
	}

	/** Stores validated arguments and shared sentinels. */
	private CalcTrueQualityParserArgsTest(final boolean fixed_, final TestFiles files_,
			final SavedState initial_){
		fixed=fixed_;
		files=files_;
		initial=initial_;
	}

	/** Runs the equals-bearing argument case and all variant-file routing controls. */
	private void process(){
		runEqualsCase();
		runVcfAliasCases();
		runPlainTextVarFileControl();
		runExplicitVcfControl();
		System.out.println("CALCTRUEQUALITY_PARSER_ARGS_PASS fixed="+fixed);
	}

	/** Verifies split-on-first semantics for input and qhist values containing '='. */
	private void runEqualsCase(){
		final Observed observed=construct("equals", files.fullInput, "qhist="+
				files.fullQhist);
		observed.checkInput("equals", fixed ? files.fullInput : files.truncatedInput);
		observed.checkQhist("equals", fixed ? files.fullQhist : files.truncatedQhist);
		observed.checkVarFile("equals", null);
		observed.checkVcfFile("equals", null);
		observed.checkReadStats("equals", true);
		System.out.println("PARSER_EQUALS_CASE_PASS fixed="+fixed);
	}

	/** Verifies every known-var alias routes a .vcf path according to phase. */
	private void runVcfAliasCases(){
		for(final String alias : VARFILE_ALIASES){
			final Observed observed=construct(alias, files.plainInput, alias+"="+
					files.vcf);
			observed.checkInput(alias, files.plainInput);
			observed.checkQhist(alias, null);
			observed.checkVarFile(alias, fixed ? null : files.vcf);
			observed.checkVcfFile(alias, fixed ? files.vcf : null);
			observed.checkReadStats(alias, false);
			System.out.println("PARSER_VARFILE_VCF_ALIAS_PASS alias="+alias+
					" fixed="+fixed);
		}
	}

	/** Verifies a non-VCF varfile remains a text-var path in both phases. */
	private void runPlainTextVarFileControl(){
		final Observed observed=construct("txt", files.plainInput,
				"varfile="+files.txt);
		observed.checkInput("txt", files.plainInput);
		observed.checkQhist("txt", null);
		observed.checkVarFile("txt", files.txt);
		observed.checkVcfFile("txt", null);
		observed.checkReadStats("txt", false);
		System.out.println("PARSER_VARFILE_TXT_CONTROL_PASS fixed="+fixed);
	}

	/** Verifies the explicit vcf option still targets the VCF field. */
	private void runExplicitVcfControl(){
		final Observed observed=construct("vcf", files.plainInput,
				"vcf="+files.vcf);
		observed.checkInput("vcf", files.plainInput);
		observed.checkQhist("vcf", null);
		observed.checkVarFile("vcf", null);
		observed.checkVcfFile("vcf", files.vcf);
		observed.checkReadStats("vcf", false);
		System.out.println("PARSER_VCF_CONTROL_PASS fixed="+fixed);
	}

	/** Constructs one fresh CalcTrueQuality and snapshots only parsed constructor fields. */
	private Observed construct(final String label, final Path input,
			final String arg){
		initial.restoreAll();
		try{
			final ArrayList<String> args=new ArrayList<String>();
			args.add("in="+input);
			args.add("bf1=t");
			args.add(arg);
			final CalcTrueQuality ctq=new CalcTrueQuality(args.toArray(new String[args.size()]));
			return new Observed(ctq);
		}catch(RuntimeException e){
			throw fail(label+" constructor failed: "+e.getMessage(), e);
		}finally{
			initial.restoreAll();
		}
	}

	/** Snapshot of the private fields parsed by the constructor. */
	private static final class Observed{
		Observed(final CalcTrueQuality ctq){
			in=strings(ctq, F_IN);
			qhist=string(ctq, F_QHIST);
			varFile=string(ctq, F_VARFILE);
			vcfFile=string(ctq, F_VCFFILE);
			readstats=object(ctq, F_READSTATS);
		}

		/** Checks the exact single input path. */
		void checkInput(final String label, final Path expected){
			check(in!=null && in.length==1, label+" in must have one entry");
			check(expected.toString().equals(in[0]), label+" in expected="+
				expected+" observed="+in[0]);
		}

		/** Checks the exact qhist path. */
		void checkQhist(final String label, final Path expected){
			checkString(label, "qhist", expected, qhist);
		}

		/** Checks the exact varFile path. */
		void checkVarFile(final String label, final Path expected){
			checkString(label, "varFile", expected, varFile);
		}

		/** Checks the exact vcfFile path. */
		void checkVcfFile(final String label, final Path expected){
			checkString(label, "vcfFile", expected, vcfFile);
		}

		/** Checks whether qhist triggered ReadStats bookkeeping allocation. */
		void checkReadStats(final String label, final boolean expected){
			check((readstats!=null)==expected, label+" readstats expected="+
				expected+" observed="+(readstats!=null));
		}

		private final String[] in;
		private final String qhist;
		private final String varFile;
		private final String vcfFile;
		private final Object readstats;
	}

	/** Readable sentinels that let the constructor validate inputs without processing. */
	private static final class TestFiles{
		/** Creates the SAM and variant sentinels. */
		static TestFiles create(){
			TestFiles result=null;
			try{
				final Path dir=Files.createTempDirectory("ctq-parser-args-");
				result=new TestFiles(dir);
				result.writeSentinels();
				return result;
			}catch(IOException e){
				if(result!=null){
					try{
						result.cleanup();
					}catch(RuntimeException cleanupFailure){
						e.addSuppressed(cleanupFailure);
					}
				}
				throw fail("Could not create parser sentinel files", e);
			}
		}

		/** Stores all sentinel paths. */
		private TestFiles(final Path dir_){
			dir=dir_;
			truncatedInput=dir.resolve("input");
			fullInput=dir.resolve("input=full.sam");
			plainInput=dir.resolve("plain.sam");
			truncatedQhist=dir.resolve("quality");
			fullQhist=dir.resolve("quality=hist.txt");
			vcf=dir.resolve("vars.vcf");
			txt=dir.resolve("vars.txt");
		}

		/** Writes readable input files used only by Tools.testInputFiles. */
		private void writeSentinels() throws IOException{
			write(truncatedInput, SAM_SENTINEL);
			write(fullInput, SAM_SENTINEL);
			write(plainInput, SAM_SENTINEL);
			write(vcf, VCF_SENTINEL);
			write(txt, TEXT_VARIANTS_SENTINEL);
		}

		/** Removes the finite sentinel tree. */
		void cleanup(){
			delete(truncatedInput);
			delete(fullInput);
			delete(plainInput);
			delete(truncatedQhist);
			delete(fullQhist);
			delete(vcf);
			delete(txt);
			delete(dir);
		}

		final Path dir;
		final Path truncatedInput;
		final Path fullInput;
		final Path plainInput;
		final Path truncatedQhist;
		final Path fullQhist;
		final Path vcf;
		final Path txt;
	}

	/** Process-global constructor and ReadStats bookkeeping restored around every case. */
	private static final class SavedState{
		SavedState(){
			parse8=SamLine.PARSE_8;
			parseOptionalMdOnly=SamLine.PARSE_OPTIONAL_MD_ONLY;
			usePigz=ReadWrite.USE_PIGZ;
			useUnpigz=ReadWrite.USE_UNPIGZ;
			zipLevel=ReadWrite.ZIPLEVEL;
			bf1=ByteFile.FORCE_MODE_BF1;
			bf2=ByteFile.FORCE_MODE_BF2;
			bf3=ByteFile.FORCE_MODE_BF3;
			bf4=ByteFile.FORCE_MODE_BF4;
			readStatsObjects=new ArrayList<ReadStats>(ReadStats.objectList);
			merged=ReadStats.merged;
		}

		/** Restores constructor statics and removes any qhist-created ReadStats objects. */
		void restoreAll(){
			SamLine.PARSE_8=parse8;
			SamLine.PARSE_OPTIONAL_MD_ONLY=parseOptionalMdOnly;
			ReadWrite.USE_PIGZ=usePigz;
			ReadWrite.USE_UNPIGZ=useUnpigz;
			ReadWrite.ZIPLEVEL=zipLevel;
			ByteFile.FORCE_MODE_BF1=bf1;
			ByteFile.FORCE_MODE_BF2=bf2;
			ByteFile.FORCE_MODE_BF3=bf3;
			ByteFile.FORCE_MODE_BF4=bf4;
			ReadStats.objectList.clear();
			ReadStats.objectList.addAll(readStatsObjects);
			ReadStats.merged=merged;
		}

		final boolean parse8;
		final boolean parseOptionalMdOnly;
		final boolean usePigz;
		final boolean useUnpigz;
		final int zipLevel;
		final boolean bf1;
		final boolean bf2;
		final boolean bf3;
		final boolean bf4;
		final ArrayList<ReadStats> readStatsObjects;
		final ReadStats merged;
	}

	/** Parses and stores CalcTrueQualityParserArgsTest arguments. */
	private static final class Args{
		Args(final String[] args0){
			final String[] args=(args0==null ? new String[0] : args0);
			final Parser parser=new Parser();
			for(final String arg : args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String key=split[0].toLowerCase(Locale.ROOT);
				final String value=(split.length>1 ? split[1] : null);
				if("fixed".equals(key)){
					fixed=Parse.parseBoolean(value);
				}else if("phase".equals(key)){
					fixed="fixed".equalsIgnoreCase(value);
					check(fixed || "baseline".equalsIgnoreCase(value),
							"phase must be baseline or fixed: "+value);
				}else if(parser.parse(arg, key, value)){
					// Consumed by BBTools common parsing.
				}else{
					throw fail("Unknown argument: "+key);
				}
			}
		}

		/** Whether parser expectations should require the repaired product code. */
		boolean fixed;
	}

	/** Creates a private-field handle. */
	private static Field field(final String name){
		try{
			final Field result=CalcTrueQuality.class.getDeclaredField(name);
			result.setAccessible(true);
			return result;
		}catch(NoSuchFieldException e){
			throw fail("Missing CalcTrueQuality field: "+name, e);
		}
	}

	/** Reads one String field. */
	private static String string(final CalcTrueQuality ctq, final Field field){
		return (String)object(ctq, field);
	}

	/** Reads one String-array field. */
	private static String[] strings(final CalcTrueQuality ctq, final Field field){
		return (String[])object(ctq, field);
	}

	/** Reads one reflected field and reports setup failures as fixture failures. */
	private static Object object(final CalcTrueQuality ctq, final Field field){
		try{return field.get(ctq);}
		catch(IllegalAccessException e){throw fail("Could not read "+field.getName(), e);}
	}

	/** Writes ASCII fixture data. */
	private static void write(final Path path, final String data) throws IOException{
		Files.write(path, data.getBytes(StandardCharsets.US_ASCII));
	}

	/** Deletes one sentinel path. */
	private static void delete(final Path path){
		try{Files.deleteIfExists(path);}
		catch(IOException e){throw fail("Could not delete "+path, e);}
	}

	/** Checks a nullable path field. */
	private static void checkString(final String label, final String field,
			final Path expected, final String observed){
		final String expectedString=(expected==null ? null : expected.toString());
		check(expectedString==null ? observed==null : expectedString.equals(observed),
				label+" "+field+" expected="+expectedString+" observed="+observed);
	}

	/** Throws when a condition is false. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw fail(message);}
	}

	/** Returns a RuntimeException with a diagnostic message. */
	private static RuntimeException fail(final String message){
		return new RuntimeException(message);
	}

	/** Returns a RuntimeException with a diagnostic message and cause. */
	private static RuntimeException fail(final String message, final Throwable cause){
		return new RuntimeException(message, cause);
	}

	private static final String[] VARFILE_ALIASES={"vars", "variants", "variations",
		"varfile", "inv"};

	private static final String SAM_SENTINEL="@HD\tVN:1.6\tSO:unsorted\n"+
		"@SQ\tSN:chr1\tLN:12\n"+
		"r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\tMD:Z:4\n";
	private static final String VCF_SENTINEL="##fileformat=VCFv4.2\n"+
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
	private static final String TEXT_VARIANTS_SENTINEL="# text variant sentinel\n";

	private static final Field F_IN=field("in");
	private static final Field F_QHIST=field("qhist");
	private static final Field F_VARFILE=field("varFile");
	private static final Field F_VCFFILE=field("vcfFile");
	private static final Field F_READSTATS=field("readstats");

	private final boolean fixed;
	private final TestFiles files;
	private final SavedState initial;
}
