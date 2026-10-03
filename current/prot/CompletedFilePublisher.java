package prot;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.LinkOption;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;

/** Publishes a closed, caller-owned regular temporary file; both paths must be on the same filesystem. */
public final class CompletedFilePublisher {

	private CompletedFilePublisher(){}

	/**
	 * With overwrite=false, createLink atomically creates the destination only if absent.
	 * ATOMIC_MOVE alone cannot provide this guarantee: its treatment of an existing target
	 * is provider-specific. Hard-link failure is propagated without a copying fallback that
	 * could expose incomplete destination bytes. The caller owns cleanup of a surviving temp.
	 *
	 * If removal of the temp fails after a successful link, the complete destination remains
	 * published and the IOException is propagated; neither path is silently deleted on error.
	 * With overwrite=true the caller explicitly authorizes atomic replacement.
	 */
	public static void publish(final Path completedTemp, final Path destination, final boolean overwrite) throws IOException {
		if(!Files.isRegularFile(completedTemp, LinkOption.NOFOLLOW_LINKS)){
			throw new IOException("Completed temporary is not a regular file: "+completedTemp);
		}
		if(completedTemp.toAbsolutePath().normalize().equals(destination.toAbsolutePath().normalize())
				|| (Files.exists(destination, LinkOption.NOFOLLOW_LINKS) && Files.isSameFile(completedTemp, destination))){
			throw new IOException("Temporary and destination name the same file: "+destination);
		}
		if(overwrite){
			Files.move(completedTemp, destination, StandardCopyOption.ATOMIC_MOVE, StandardCopyOption.REPLACE_EXISTING);
		}else{
			Files.createLink(destination, completedTemp);
			Files.delete(completedTemp);
		}
	}
}
