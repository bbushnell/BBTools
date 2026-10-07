package stream;

import java.util.ArrayList;

import fileIO.FileFormat;

/**
 * File-backed target for MultiFileWriter using the explicit lightweight factory.
 * Retains descriptors and header references, not open handles. Each reopening
 * preserves explicit format/compression overrides and switches to append.
 * Caller must map distinct logical destinations to disjoint physical outputs and
 * keep header/dictionary arrays unchanged for the target's lifetime.
 *
 * @author Shinobu
 * @date October 1, 2026
 */
public final class FileWriterTarget implements MultiFileWriter.Target{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates a target without opening files.
	 * @param factory_ Shared immutable leaf factory
	 * @param first_ Primary output descriptor
	 * @param second_ Optional separate mate output descriptor
	 * @param qual1_ Optional primary FASTA quality filename
	 * @param qual2_ Optional mate FASTA quality filename
	 * @param header_ Optional immutable SAM/BAM header reference
	 * @param useSharedHeader_ Use the shared input header instead of header_ */
	public FileWriterTarget(final LightweightWriterFactory factory_, final FileFormat first_, final FileFormat second_,
			final String qual1_, final String qual2_, final ArrayList<byte[]> header_, final boolean useSharedHeader_){
		if(factory_==null || first_==null || (second_==null && qual2_!=null)){
			throw new IllegalArgumentException("A factory and primary output are required; mate QUAL needs a mate output");
		}
		factory=factory_;
		first=first_;
		second=second_;
		qual1=qual1_;
		qual2=qual2_;
		header=header_;
		useSharedHeader=useSharedHeader_;
		outputs=1+(second==null ? 0 : 1)+(qual1==null ? 0 : 1)+(qual2==null ? 0 : 1);
		reopenable=reopenable(first) && reopenable(second)
			&& reopenable(qual1==null ? null : LightweightWriterFactory.qualFormat(first, qual1))
			&& reopenable(qual2==null ? null : LightweightWriterFactory.qualFormat(second, qual2));
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** @return Physical outputs, including separate mates and quality files */
	@Override
	public int outputCount(){return outputs;}

	/** @return False for ZIP or non-file targets, which cannot rotate safely here */
	@Override
	public boolean canReopen(){return reopenable;}

	/** Opens a new leaf generation with stable serialization and dense local IDs.
	 * @param reopen Whether an earlier generation was closed
	 * @return Unstarted ordinary Writer */
	@Override
	public Writer open(final boolean reopen){
		if(reopen && !reopenable){throw new IllegalStateException("Target cannot be reopened for append: "+first.name());}
		return factory.getStream(descriptor(first, reopen), descriptor(second, reopen), qual1, qual2, header, useSharedHeader);
	}

	/** Preserves descriptor overrides; ST uses ordered local IDs, ZT needs no ordering queue. */
	private FileFormat descriptor(final FileFormat ff, final boolean reopen){
		return ff==null ? null : FileFormat.testOutput(ff.name(), ff.format(), ff.format(), ff.compression(), false,
			ff.overwrite(), ff.append() || reopen, factory.mode==MultiWriterPolicy.ST);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates a resolver for legacy-style percent-substitution output patterns.
	 * Every nonnull pattern must contain a percent sign; all occurrences are replaced
	 * literally by the destination name. Each resulting filename selects its own
	 * format. Callers must choose patterns/names producing disjoint physical files.
	 * @param factory Shared lightweight leaf factory
	 * @param first Primary output pattern
	 * @param second Optional mate output pattern
	 * @param qual1 Optional primary QUAL pattern
	 * @param qual2 Optional mate QUAL pattern
	 * @param defaultFormat Format fallback for filenames without a recognized extension
	 * @param overwrite Permit replacement on the first open
	 * @param append Append on the first open as well as subsequent generations
	 * @param header Borrowed immutable SAM/BAM header
	 * @param sharedHeader Use the shared input header instead of header
	 * @return Resolver that opens no files until the container submits output
	 */
	public static MultiFileWriter.TargetFactory patternFactory(final LightweightWriterFactory factory,
			final String first, final String second, final String qual1, final String qual2,
			final int defaultFormat, final boolean overwrite, final boolean append,
			final ArrayList<byte[]> header, final boolean sharedHeader){
		if(factory==null || first==null || (second==null && qual2!=null)){
			throw new IllegalArgumentException("Patterns require a factory and primary output; mate QUAL needs a mate output");
		}
		checkPattern(first); checkPattern(second); checkPattern(qual1); checkPattern(qual2);
		return name->new FileWriterTarget(factory,
			FileFormat.testOutput(expand(first, name), defaultFormat, null, false, overwrite, append, false),
			FileFormat.testOutput(expand(second, name), defaultFormat, null, false, overwrite, append, false),
			expand(qual1, name), expand(qual2, name), header, sharedHeader);
	}

	/** Requires substitution so different logical names do not trivially overwrite one file. */
	private static void checkPattern(final String pattern){
		if(pattern!=null && pattern.indexOf('%')<0){throw new IllegalArgumentException("Output pattern needs a percent placeholder: "+pattern);}
	}

	/** Expands a nullable pattern without regular-expression replacement semantics. */
	private static String expand(final String pattern, final String name){return pattern==null ? null : pattern.replace("%", name);}

	/** Determines whether the supported compression/container can append another generation. */
	private static boolean reopenable(final FileFormat ff){
		return ff==null || (ff.file() && ff.compression()!=FileFormat.ZIP);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Shared immutable leaf execution policy. */
	private final LightweightWriterFactory factory;
	/** Original descriptors, including explicit format/compression overrides. */
	private final FileFormat first, second;
	/** Optional numeric quality output names. */
	private final String qual1, qual2;
	/** Borrowed immutable dictionary/header, reused across generations. */
	private final ArrayList<byte[]> header;
	/** Header source choice and whether all physical outputs support reopening. */
	private final boolean useSharedHeader, reopenable;
	/** Physical outputs per open leaf generation. */
	private final int outputs;
}
