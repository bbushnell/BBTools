package prok;

/** Optional generated-window observer. A group is not a scheduled alignment.
 * Callbacks may retain or mutate their owned member arrays, never caller state.
 * Shared implementations must serialize concurrent caller threads.
 * @author Keqing
 */
interface NcrnaVoteDiagSink {
	void group(NcrnaVoteDiagnostics.Group group);
}
