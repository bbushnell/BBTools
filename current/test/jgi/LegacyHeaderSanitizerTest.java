package test.jgi;

import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;

import jgi.ReformatReads;
import jgi.ReformatReadsOld;
import processor.ReformatProcessor;

/**
 * Literal direct String-sanitizer controls for the modern and two legacy helpers.
 * ASCII mapping follows the established modern table. Every non-ASCII UTF-16
 * unit maps to one underscore; a supplementary pair therefore maps to two.
 * Baseline legacy calls instead expose their out-of-table access. Null, empty
 * and unchanged headers retain their input reference. Run once in a fresh JVM
 * with UTF-8 selected explicitly; no Read, driver, parser or worker is created.
 * @author Sucrose
 */
public final class LegacyHeaderSanitizerTest{

	/** Requires phase=baseline|fixed; the modern counterpart must pass in either phase. */
	public static void main(final String[] args){
		if(args.length!=1 || (!args[0].equals("phase=baseline") && !args[0].equals("phase=fixed"))){
			throw new IllegalArgumentException("Set phase=baseline|fixed.");
		}
		final boolean baseline=args[0].equals("phase=baseline");
		final Charset charset=Charset.defaultCharset();
		check(charset.equals(StandardCharsets.UTF_8),
			"Sanitizers construct Strings from ASCII bytes using the default charset; this fixture requires UTF-8: "+charset);
		System.out.println("header_sanitizer_encoding\tdefault="+charset.name()+
			"\tfile.encoding="+System.getProperty("file.encoding"));
		int failures=0, references=0, cases=0;
		for(int variant=0; variant<VARIANTS.length; variant++){
			for(int i=0; i<CASES.length; i++){
				final String[] row=CASES[i];
				final String label=VARIANTS[variant]+"-"+row[0];
				final String input=row[1]==null ? null : new String(row[1]);
				final boolean expectFailure=baseline && variant>0 && i>=ASCII_CASES;
				String actual=null;
				ArrayIndexOutOfBoundsException failure=null;
				try{actual=convert(variant, input);}
				catch(ArrayIndexOutOfBoundsException e){failure=e;}
				check((failure!=null)==expectFailure,
					label+": expected array-index failure="+expectFailure+", observed="+(failure!=null));
				if(failure!=null){
					final StackTraceElement top=failure.getStackTrace()[0];
					final String owner=variant==1 ? "jgi.ReformatReads" : "jgi.ReformatReadsOld";
					check(top.getClassName().equals(owner) && top.getMethodName().equals("fixHeader"),
						label+": failure must originate in the selected legacy sanitizer: "+top);
					failures++;
				}else{
					final boolean same=actual==input;
					final boolean equal=row[2]==null ? actual==null : row[2].equals(actual);
					check(equal && same==(i<IDENTITY_CASES),
						label+": expected "+display(row[2])+", observed "+display(actual)+", same reference="+same);
					if(input!=null && same){references++;}
				}
				System.out.println("header_sanitizer_case\t"+label+"\tinput="+display(input)+
					"\tintended="+display(row[2])+"\toutcome="+(failure==null ? "value" : "array-index-failure")+
					"\tactual="+(failure==null ? display(actual) : "<absent>")+
					"\tsameReference="+(failure==null ? Boolean.toString(actual==input) : "NA"));
				cases++;
			}
		}
		check(failures==(baseline ? 16 : 0),
			"Eight non-ASCII cases per legacy helper must fail only before the guard: "+failures);
		check(references==6, "Empty and allowed-only nonnull Strings must retain identity for all three helpers: "+references);
		if(checks!=99 || cases!=48){
			throw new AssertionError("Frozen sanitizer fixture requires 99 checks/48 cases: "+checks+"/"+cases);
		}
		System.out.println("LEGACY_HEADER_SANITIZER_PASS phase="+(baseline ? "baseline" : "fixed")+
			" checks="+checks+" cases="+cases+" arrayIndexFailures="+failures+" unchangedReferences="+references+" modernControls=16");
	}

	/** Uses actual public static String overloads; this dispatch is only functional setup. */
	private static String convert(final int variant, final String header){
		assert(variant>=0 && variant<VARIANTS.length) :
			"The selected helper and diagnostic label must use the same three-entry variant table: "+variant;
		if(variant==0){return ReformatProcessor.fixHeader(header);}
		if(variant==1){return ReformatReads.fixHeader(header);}
		return ReformatReadsOld.fixHeader(header);
	}

	/** Covers every ASCII table index without deriving expected mappings from either helper. */
	private static String allAscii(){
		final char[] chars=new char[128];
		for(int i=0; i<chars.length; i++){chars[i]=(char)i;}
		return new String(chars);
	}

	/** ASCII-only observations retain exact UTF-16 units and avoid tab/newline ambiguity. */
	private static String display(final String value){
		if(value==null){return "<null>";}
		final StringBuilder out=new StringBuilder("[");
		for(int i=0; i<value.length(); i++){
			final char c=value.charAt(i);
			if(c>=32 && c<127 && c!='\\'){out.append(c);}
			else{
				out.append("\\u");
				for(int shift=12; shift>=0; shift-=4){out.append(HEX[(c>>>shift)&15]);}
			}
		}
		return out.append(']').toString();
	}

	/** Literal checks stay active even if JVM assertions are accidentally disabled. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
		checks++;
	}

	/** Label, input, intended output. First eight are ASCII/null; first three retain identity. */
	private static final String[][] CASES={
		{"null", null, null},
		{"empty", "", ""},
		{"allowed", "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_.#-()~",
			"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_.#-()~"},
		{"space", "A B", "A_B"},
		{"separators", "a/b\\c:d", "a_b_c_d"},
		{"controls", "\000\t\n\r", "____"},
		{"ascii-127", "a\u007fb", "a_b"},
		{"all-ascii", allAscii(),
			"________________"+"________________"+"___#____()___-._"+"0123456789______"+
			"_ABCDEFGHIJKLMNO"+"PQRSTUVWXYZ_____"+"_abcdefghijklmno"+"pqrstuvwxyz___~_"},
		{"boundary-128", "a\u0080b", "a_b"},
		{"latin", "caf\u00e9", "caf_"},
		{"cjk", "A\u6f22B", "A_B"},
		{"supplementary", "A\ud83d\ude00B", "A__B"},
		{"multiple-nonascii", "\u0080\u00e9\u6f22", "___"},
		{"allowed-and-nonascii", "._#-()~\u0080AZ09", "._#-()~_AZ09"},
		{"boundary-127-128", "\u007f\u0080", "__"},
		{"utf16-maximum", "\uffff", "_"}
	};
	private static final String[] VARIANTS={"processor", "legacy", "legacyold"};
	private static final int ASCII_CASES=8, IDENTITY_CASES=3;
	private static final char[] HEX="0123456789abcdef".toCharArray();
	private static int checks;
}
