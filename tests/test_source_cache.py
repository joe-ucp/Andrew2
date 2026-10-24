"""Verify remote-only source retrieval, integrity checks and offline propagation."""
import hashlib
import io
import os
from pathlib import Path
import runpy
import sys
import tempfile
import unittest
from unittest import mock
import urllib.error

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'code'))
import fetch_inputs


class SourceCacheTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='beta-cell-cache-test-')
        self.addCleanup(self.temporary.cleanup)
        self.cache = Path(self.temporary.name) / 'cache'
        self.payload = b'Synthetic source-cache fixture; no third-party code.'
        patcher = mock.patch.dict(fetch_inputs.ROTTNER_ARCHIVE,
                                  {'sha256': hashlib.sha256(self.payload).hexdigest()})
        patcher.start()
        self.addCleanup(patcher.stop)

    def populate(self, directory=None, payload=None):
        directory = directory or self.cache
        directory.mkdir(parents=True, exist_ok=True)
        path = directory / fetch_inputs.ROTTNER_ARCHIVE['file']
        path.write_bytes(self.payload if payload is None else payload)
        return path

    def test_download_then_offline_reuse(self):
        with mock.patch.object(fetch_inputs.urllib.request, 'urlopen',
                               return_value=io.BytesIO(self.payload)) as download:
            path = fetch_inputs.resolve_rottner_archive(self.cache)
            self.assertEqual(path.read_bytes(), self.payload)
            self.assertEqual(download.call_args.args[0].full_url,
                             fetch_inputs.ROTTNER_ARCHIVE['url'])
            self.assertEqual(list(self.cache.glob('*.part')), [])
        with mock.patch.object(fetch_inputs.urllib.request, 'urlopen') as download:
            self.assertEqual(fetch_inputs.resolve_rottner_archive(self.cache, offline=True), path)
            download.assert_not_called()

    def test_offline_missing_is_actionable_and_never_downloads(self):
        with mock.patch.object(fetch_inputs.urllib.request, 'urlopen') as download:
            with self.assertRaisesRegex(FileNotFoundError, 'fetch_inputs.py --source-cache'):
                fetch_inputs.resolve_rottner_archive(self.cache, offline=True)
            download.assert_not_called()
        self.assertFalse(self.cache.exists())

    def test_corrupt_cache_is_rejected_without_overwriting(self):
        path = self.populate(payload=b'corrupt')
        for offline in (True, False):
            with self.subTest(offline=offline), mock.patch.object(fetch_inputs.urllib.request, 'urlopen') as download:
                with self.assertRaisesRegex(ValueError, 'HASH MISMATCH for cached'):
                    fetch_inputs.resolve_rottner_archive(self.cache, offline=offline)
                download.assert_not_called()
                self.assertEqual(path.read_bytes(), b'corrupt')

    def test_corrupt_download_is_not_cached(self):
        with mock.patch.object(fetch_inputs.urllib.request, 'urlopen', return_value=io.BytesIO(b'corrupt')):
            with self.assertRaisesRegex(ValueError, 'nothing was cached'):
                fetch_inputs.resolve_rottner_archive(self.cache)
        self.assertFalse(self.cache.exists())

    def test_network_failure_is_actionable(self):
        with mock.patch.object(fetch_inputs.urllib.request, 'urlopen', side_effect=urllib.error.URLError('unavailable')):
            with self.assertRaisesRegex(RuntimeError, 'verified external cache via --source-cache'):
                fetch_inputs.resolve_rottner_archive(self.cache)
        self.assertFalse(self.cache.exists())

    def test_repository_cache_is_rejected_before_network(self):
        with mock.patch.object(fetch_inputs.urllib.request, 'urlopen') as download:
            with self.assertRaisesRegex(ValueError, 'outside the repository'):
                fetch_inputs.resolve_rottner_archive(ROOT / 'source-cache')
            download.assert_not_called()

    def test_explicit_cache_overrides_environment(self):
        explicit = self.populate()
        other = self.populate(Path(self.temporary.name) / 'environment-cache')
        with mock.patch.dict(os.environ, {'BETA_CELL_SOURCE_CACHE': str(other.parent)}):
            self.assertEqual(fetch_inputs.resolve_rottner_archive(offline=True), other.resolve())
            self.assertEqual(fetch_inputs.resolve_rottner_archive(self.cache, offline=True), explicit.resolve())

    def test_default_cache_is_under_user_home(self):
        home = Path(self.temporary.name) / 'home'
        path = self.populate(home / '.cache' / 'beta-cell-screen-comparison' / 'sources')
        with mock.patch.dict(os.environ, {'BETA_CELL_SOURCE_CACHE': ''}), mock.patch.object(Path, 'home', return_value=home):
            self.assertEqual(fetch_inputs.resolve_rottner_archive(offline=True), path.resolve())

    def test_reproduction_propagates_cache_and_offline_options(self):
        args = ['reproduce.py', '--offline', '--source-cache', str(self.cache), '--build-manuscript']
        with mock.patch.object(sys, 'argv', args), mock.patch('subprocess.run') as run, mock.patch('builtins.print'):
            runpy.run_path(str(ROOT / 'reproduce.py'), run_name='__main__')
        commands = [call.args[0][1:] for call in run.call_args_list]
        source_commands = [c for c in commands if c[0] in ('code/fetch_inputs.py', 'code/audit_source_semantics.py')]
        self.assertEqual(len(source_commands), 3)
        for command in source_commands:
            self.assertIn('--offline', command)
            self.assertEqual(command[command.index('--source-cache') + 1], str(self.cache))
        validation = next(c for c in commands if c[0] == 'tests/validate_source_semantics.py')
        self.assertEqual(validation[validation.index('--source-cache') + 1], str(self.cache))
        self.assertEqual(commands[-1], ['code/build_manuscript.py', '--compile'])


if __name__ == '__main__':
    unittest.main()
