#include "targetver.h"

#include <assert.h>
#include <stdlib.h>
#include <windows.h>
#include "resource.h"
#include <commdlg.h>
#include <shlobj.h>

#define ARRAY_SZ(arr) (sizeof(arr) / sizeof(arr[0]))
#define WLEN(l) ((l) * sizeof(WCHAR))

#define HAlloc(sz, h) HeapAlloc(h, HEAP_GENERATE_EXCEPTIONS, sz)
#define HReAlloc(ptr, sz, h) HeapReAlloc(h, HEAP_GENERATE_EXCEPTIONS, ptr, sz)
#define HFree(ptr, h) HeapFree(h, 0, ptr)

#define UI_TOP_OFFSET 16
#define UI_CAPTION_SEP 8
#define UI_BTN_HEIGHT 32
#define UI_BTN_WIDTH 128
#define UI_HGT_PAD 8
#define UI_SIDE_OFFSET 8

#define MMB_RUNTIME_DEFAULT_PATH L"MMB_rt"
#define MMB_MODULE L"MMB"
#define MMB_EXEC L"MMB.exe"
#define MMB_PARAMS_FILE L"parameters.csv"

#define WM_MMB_QUIT (WM_USER + 1)
#define WM_MMB_OUTPUT (WM_USER + 2)

#define MMB_ST_NOT_STARTED_CAP L"Not started"
#define MMB_ST_RUNNING_CAP L"Running"
#define MMB_ST_FAILED_CAP L"Failed"
#define MMB_ST_FINISHED_CAP L"Finished"

#define CloseInvalidateHandle(h) { CloseHandle(h); h = NULL; }

enum MmbStatus {
	MMB_ST_NOT_STARTED,
	MMB_ST_RUNNING,
	MMB_ST_FAILED,
	MMB_ST_FINISHED
};

static const WCHAR g_szWindowTitle[] = L"WinMMB";
static const WCHAR g_szWindowClass[] = L"WinMMBWndClass";

static const WCHAR g_szMmbRtDirPathCap[] = L"MMB runtime path";
static const WCHAR g_szMmbCmdsPathCap[] = L"MMB commands file";
static const WCHAR g_szMmbWorkDirCap[] = L"MMB working directory";
static const WCHAR g_szMmbOutputCap[] = L"MMB output";
static const WCHAR g_szMmbStatusCap[] = L"MMB status";

static const LPCWSTR g_formCaps[] = {
	g_szMmbRtDirPathCap,
	g_szMmbCmdsPathCap,
	g_szMmbWorkDirCap,
	g_szMmbStatusCap
};

static LONG g_captionWidest = 0;
static LONG g_captionHighest = 0;
static LONG g_formStep = 0;

static LPWSTR g_mmbRtPath = NULL;
static LPWSTR g_cmdsPath = NULL;
static LPWSTR g_mmbWorkDir = NULL;
static LPWSTR g_mmbOutput = NULL;
static BOOL   g_scrollMmbOutput = TRUE;

static HFONT g_hMainFont = NULL;
static HWND g_hMmbRtDirCap = NULL;
static HWND g_hMmbCmdsPathCap = NULL;
static HWND g_hMmbWorkDirCap = NULL;
static HWND g_hMmbOutputCap = NULL;
static HWND g_hMmbRtDir = NULL;
static HWND g_hMmbCmdsPath = NULL;
static HWND g_hMmbWorkDir = NULL;
static HWND g_hMmbOutput = NULL;
static HWND g_hMmbRtDirBrowse = NULL;
static HWND g_hMmbCmdsPathBrowse = NULL;
static HWND g_hOutputDirBrowse = NULL;
static HWND g_hRunBtn = NULL;
static HWND g_hStopBtn = NULL;
static HWND g_hMmbStatusCap = NULL;
static HWND g_hMmbStatus = NULL;
static HWND g_hDontScroll = NULL;
static HWND g_hMmbOutputToFile = NULL;

static PROCESS_INFORMATION g_mmbProc;
static HANDLE g_hMmbReadPipe = NULL;
static HANDLE g_hMonitorThread = NULL;

static HANDLE g_hOutputMutex = NULL;

static LPVOID g_heap = NULL;

enum {
	ID_MMB_RT = 1,
	ID_MMB_CMDS_PATH,
	ID_OUTPUT_DIR,
	ID_MMB_OUTPUT,

	IDC_BROWSE_MMB_RT,
	IDC_BROWSE_CMDS_PATH,
	IDC_BROWSE_OUTPUT_DIR,
	IDC_RUN_MMB,
	IDC_STOP_MMB,
	IDC_DONT_SCROLL,
	IDC_MMB_OUTPUT_TO_FILE
};

static
VOID FreeAllocated()
{
	HFree(g_mmbOutput, g_heap);
	HFree(g_mmbRtPath, g_heap);
	HFree(g_cmdsPath, g_heap);
	HFree(g_mmbWorkDir, g_heap);
}

static
BOOL InitMainFont()
{
	NONCLIENTMETRICSW ncm = { sizeof(NONCLIENTMETRICSW) };
	if (!SystemParametersInfoW(SPI_GETNONCLIENTMETRICS, ncm.cbSize, &ncm, 0)) {
		MessageBoxW(NULL, L"Failed to get system GUI font", L"Critical error", MB_OK | MB_ICONERROR);
		return FALSE;
	}

	g_hMainFont = CreateFontIndirectW(&(ncm.lfMessageFont));
	return TRUE;
}

static
LPWSTR PrettyErrorString(LPCWSTR prefix, DWORD error)
{
	LPWSTR buf = NULL;
	LPWSTR completeMessage = NULL;

	if (FormatMessageW(
		FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
		NULL,
		error,
		MAKELANGID(LANG_NEUTRAL, SUBLANG_NEUTRAL),
		(LPWSTR) &buf,
		0,
		NULL
	)) {
		const size_t len = wcslen(prefix) + wcslen(buf) + 2;
		completeMessage = HAlloc(WLEN(len + 1), g_heap);

		wcscpy_s(completeMessage, len + 1, prefix);
		wcscat_s(completeMessage, len + 1, L": ");
		wcscat_s(completeMessage, len + 1, buf);
	} else {
		LPCWSTR unknown = L"Unknown error";

		const size_t len = wcslen(prefix) + wcslen(unknown) + 2;
		completeMessage = HAlloc(WLEN(len + 1), g_heap);

		wcscpy_s(completeMessage, len + 1, prefix);
		wcscat_s(completeMessage, len + 1, L": ");
		wcscat_s(completeMessage, len + 1, unknown);
	}

	if (buf)
		LocalFree(buf);
	return completeMessage;
}

static
LPSTR WideCharToSystemAnsi(LPCWSTR in)
{
	int count;

	count = WideCharToMultiByte(
		CP_ACP,
		WC_COMPOSITECHECK,
		in,
		-1,
		NULL,
		0,
		NULL,
		NULL
	);
	if (count < 1)
		return NULL;

	LPSTR ret = HAlloc(count + 1, g_heap);
	count = WideCharToMultiByte(
		CP_ACP,
		WC_COMPOSITECHECK,
		in,
		-1,
		ret,
		count,
		NULL,
		NULL
	);
	if (count < 1) {
		HFree(ret, g_heap);
		return NULL;
	}
	return ret;
}

static
VOID AppendMmbOutput(HWND hWnd, LPCSTR data, DWORD len)
{
	if (len < 1)
		return;

	if (WaitForSingleObject(g_hOutputMutex, INFINITE) != WAIT_OBJECT_0)
		return;

	LPWSTR buf;
        size_t cvt;
	size_t totalLen;

	buf = HAlloc(WLEN(len + 1), g_heap);

	mbstowcs_s(&cvt, buf, len + 1, data, len);

	totalLen = wcslen(g_mmbOutput) + wcslen(buf);
        if (totalLen < 1) {
                HFree(buf, g_heap);
		ReleaseMutex(g_hOutputMutex);
                return;
        }

	g_mmbOutput = HReAlloc(g_mmbOutput, WLEN(totalLen + 1), g_heap);
	wcscat_s(g_mmbOutput, totalLen + 1, buf);

	HFree(buf, g_heap);

	ReleaseMutex(g_hOutputMutex);

	PostMessageW(hWnd, WM_MMB_OUTPUT, 1, 1);
}

static
LPSTR BuildMmbCommandLine(LPWSTR cmdsPath, LPWSTR workDir)
{
	LPWSTR buf;
	DWORD len;
	DWORD lenWorkDir;

	len = wcslen(cmdsPath) + wcslen(MMB_MODULE L" -C ");
	lenWorkDir = wcslen(workDir);

	if (lenWorkDir > 0)
		len += wcslen(L" -workdir ") + lenWorkDir;

	buf = HAlloc(WLEN(len + 1), g_heap);
	wcscpy_s(buf, len + 1, MMB_MODULE L" -C ");
	wcscat_s(buf, len + 1, cmdsPath);

	if (lenWorkDir > 0) {
		wcscat_s(buf, len + 1, L" -workdir ");
		wcscat_s(buf, len + 1, workDir);
	}

	LPSTR ret = WideCharToSystemAnsi(buf);
	HFree(buf, g_heap);

	return ret;
}

static
BOOL CheckFileExists(LPCWSTR path)
{
	HANDLE hFile = CreateFileW(
		path,
		GENERIC_READ,
		0,
		NULL,
		OPEN_EXISTING,
		FILE_ATTRIBUTE_NORMAL,
		NULL
	);
	if (hFile == INVALID_HANDLE_VALUE)
		return FALSE;

	CloseHandle(hFile);

	return TRUE;
}

static
BOOL CheckMmbRtDirectory(LPCWSTR dir)
{
	const DWORD baseLen = wcslen(dir);
	LPWSTR path;
	BOOL ret;

	DWORD len = baseLen + wcslen(MMB_EXEC) + 1;
	path = HAlloc(WLEN(len + 1), g_heap);
	wcscpy_s(path, len + 1, dir);
	wcscat_s(path, len + 1, L"\\");
	wcscat_s(path, len + 1, MMB_EXEC);

	if (!CheckFileExists(path)) {
		ret = FALSE;
		goto out;
	}

	len = baseLen + wcslen(MMB_PARAMS_FILE) + 1;
	path = HReAlloc(path, WLEN(len + 1), g_heap);
	wcscpy_s(path, len + 1, dir);
	wcscat_s(path, len + 1, L"\\");
	wcscat_s(path, len + 1, MMB_PARAMS_FILE);

	if (!CheckFileExists(path)) {
		ret = FALSE;
		goto out;
	}

	ret = TRUE;
out:
	HFree(path, g_heap);

	return ret;
}

static
BOOL CopyMmbParametersFile(LPCWSTR rtDir, LPCWSTR workDir, DWORD *error)
{
	size_t len;
	LPWSTR src;
	LPWSTR dst;
	BOOL ret;

	len = wcslen(rtDir) + wcslen(MMB_PARAMS_FILE) + 1;
	src = HAlloc(WLEN(len + 1), g_heap);
	wcscpy_s(src, len + 1, rtDir);
	wcscat_s(src, len + 1, L"\\");
	wcscat_s(src, len + 1, MMB_PARAMS_FILE);

	len = wcslen(workDir) + wcslen(MMB_PARAMS_FILE) + 1;
	dst = HAlloc(WLEN(len + 1), g_heap);
	wcscpy_s(dst, len + 1, workDir);
	wcscat_s(dst, len + 1, L"\\");
	wcscat_s(dst, len + 1, MMB_PARAMS_FILE);


	HANDLE hFile = CreateFileW(
		dst,
		GENERIC_READ,
		0,
		NULL,
		OPEN_EXISTING,
		FILE_ATTRIBUTE_NORMAL,
		NULL
	);
	if (hFile != INVALID_HANDLE_VALUE) {
		CloseHandle(hFile);
		ret = TRUE;
		goto out;
	}
	CloseHandle(hFile);

	ret = CopyFileW(src, dst, FALSE);
	if (!ret)
		*error = GetLastError();

out:
	HFree(src, g_heap);
	HFree(dst, g_heap);

	return ret;
}

static
VOID ClearMmbOutput(HWND hWnd) {
	if (WaitForSingleObject(g_hOutputMutex, INFINITE) != WAIT_OBJECT_0)
		return;

	g_mmbOutput = HReAlloc(g_mmbOutput, WLEN(1), g_heap);
	g_mmbOutput[0] = 0;

	ReleaseMutex(g_hOutputMutex);

	PostMessageW(hWnd, WM_MMB_OUTPUT, 0, 0);
}

static
VOID SetMmbStatus(enum MmbStatus status)
{
	switch (status) {
	case MMB_ST_NOT_STARTED:
		SetWindowTextW(g_hMmbStatus, MMB_ST_NOT_STARTED_CAP);
		break;
	case MMB_ST_RUNNING:
		SetWindowTextW(g_hMmbStatus, MMB_ST_RUNNING_CAP);
		break;
	case MMB_ST_FAILED:
		SetWindowTextW(g_hMmbStatus, MMB_ST_FAILED_CAP);
		break;
	case MMB_ST_FINISHED:
		SetWindowTextW(g_hMmbStatus, MMB_ST_FINISHED_CAP);
		break;
	}
}

static
VOID FinalizeMmb(HWND hWnd, DWORD exitCode)
{
	CloseInvalidateHandle(g_mmbProc.hProcess);
	CloseInvalidateHandle(g_mmbProc.hThread);

	if (g_hMonitorThread == NULL)
		return;

	if (WaitForSingleObject(g_hMonitorThread, INFINITE) != WAIT_OBJECT_0) {
		MessageBoxW(
			hWnd,
			L"Failed to terminate monitor thread. This means that something "
			L"has gone very, very wrong. All I can do is commit digital suicide.",
			L"Critical error",
			MB_OK | MB_ICONERROR
		);
		abort();
	}

	CloseInvalidateHandle(g_hMonitorThread);

	SetMmbStatus(exitCode == 0 ? MMB_ST_FINISHED : MMB_ST_FAILED);
}

static
LPSTR GetMmbExecutable(LPCWSTR basePath)
{
	DWORD len = wcslen(basePath) + wcslen(MMB_EXEC) + 1;
	LPWSTR buf = HAlloc(WLEN(len + 1), g_heap);

	wcscpy_s(buf, len + 1, basePath);
	wcscat_s(buf, len + 1, L"\\");
	wcscat_s(buf, len + 1, MMB_EXEC);

	LPSTR ret = WideCharToSystemAnsi(buf);
	HFree(buf, g_heap);

	return ret;
}

static
BOOL InitializeDefaults()
{
	WCHAR buf[MAX_PATH];

	DWORD len = GetCurrentDirectoryW(MAX_PATH, buf);
	if (len < 1)
		return FALSE;

	len += wcslen(MMB_RUNTIME_DEFAULT_PATH) + 1;

	// Set default path to MMB runtime
	g_mmbRtPath = HAlloc(WLEN(len + 1), g_heap);
	g_mmbRtPath[0] = 0;

	wcscat_s(g_mmbRtPath, len + 1, buf);
	wcscat_s(g_mmbRtPath, len + 1, L"\\");
	wcscat_s(g_mmbRtPath, len + 1, MMB_RUNTIME_DEFAULT_PATH);

	// Set default path to commands file
	g_cmdsPath = HAlloc(WLEN(1), g_heap);
	g_cmdsPath[0] = 0;

	// Set default output path
	g_mmbWorkDir = HAlloc(WLEN(1), g_heap);
	g_mmbWorkDir[0] = 0;

        g_mmbOutput = HAlloc(WLEN(1), g_heap);
        g_mmbOutput[0] = 0;

	return TRUE;
}

static
LONG TextWidth(HWND hWnd, LPCWSTR text)
{
	HDC hDc = GetDC(hWnd);
	if (hDc == NULL)
		return 0;

	SelectObject(hDc, g_hMainFont);

	SIZE sz;
	GetTextExtentPointW(hDc, text, wcslen(text), &sz);
	ReleaseDC(hWnd, hDc);

	return sz.cx;
}

static
BOOL InitControls(HWND hWnd, HINSTANCE hInstance)
{
	assert(g_formStep > 0);

	RECT rect;
	if (GetClientRect(hWnd, &rect) == FALSE)
		return FALSE;

	const LONG xL = UI_SIDE_OFFSET;
	const LONG xR = xL + g_captionWidest + UI_CAPTION_SEP;
	LONG y = UI_TOP_OFFSET;

	// MMB rubtime path line
	g_hMmbRtDirCap = CreateWindowW(
		L"STATIC",
		g_szMmbRtDirPathCap,
		WS_VISIBLE | WS_CHILD | SS_CENTER,
		xL, y + (UI_HGT_PAD / 2) - 1, TextWidth(hWnd, g_szMmbRtDirPathCap), g_captionHighest,
		hWnd,
		NULL,
		hInstance,
		NULL
	);
	if (g_hMmbRtDirCap == NULL)
		return FALSE;
	SendMessageW(g_hMmbRtDirCap, WM_SETFONT, (WPARAM)g_hMainFont, TRUE);

	g_hMmbRtDir = CreateWindowW(
		L"EDIT",
		g_mmbRtPath,
		WS_BORDER | WS_CHILD | WS_VISIBLE | ES_LEFT,
		xR,  y, rect.right - xR - UI_BTN_WIDTH - UI_CAPTION_SEP - UI_SIDE_OFFSET, g_captionHighest,
		hWnd,
		(HMENU)ID_MMB_RT,
		hInstance,
		NULL
	);
	if (g_hMmbRtDir == NULL)
		return FALSE;
	SendMessageW(g_hMmbRtDir, WM_SETFONT, (WPARAM)g_hMainFont, TRUE);

	g_hMmbRtDirBrowse = CreateWindowW(
		L"BUTTON",
		L"Browse...",
		WS_TABSTOP | WS_VISIBLE | WS_CHILD | BS_DEFPUSHBUTTON,
		rect.right - UI_BTN_WIDTH - UI_SIDE_OFFSET, y, UI_BTN_WIDTH, g_captionHighest,
		hWnd,
		(HMENU)IDC_BROWSE_MMB_RT,
		hInstance,
		NULL);
	if (g_hMmbRtDirBrowse == NULL)
		return FALSE;
	SendMessageW(g_hMmbRtDirBrowse, WM_SETFONT, (WPARAM)g_hMainFont, TRUE);

	// MMB commands line
	y += g_formStep;

	g_hMmbCmdsPathCap = CreateWindowW(
		L"STATIC",
		g_szMmbCmdsPathCap,
		WS_VISIBLE | WS_CHILD | SS_CENTER,
		xL, y + (UI_HGT_PAD / 2) - 1, TextWidth(hWnd, g_szMmbCmdsPathCap), g_captionHighest,
		hWnd,
		NULL,
		hInstance,
		NULL
	);
	if (g_hMmbCmdsPathCap == NULL)
		return FALSE;
	SendMessageW(g_hMmbCmdsPathCap, WM_SETFONT, (WPARAM)g_hMainFont, TRUE);

	g_hMmbCmdsPath = CreateWindowW(
		L"EDIT",
		g_cmdsPath,
		WS_BORDER | WS_CHILD | WS_VISIBLE | ES_LEFT,
		xR,  y, rect.right - xR - UI_BTN_WIDTH - UI_CAPTION_SEP - UI_SIDE_OFFSET, g_captionHighest,
		hWnd,
		(HMENU)ID_MMB_CMDS_PATH,
		hInstance,
		NULL
	);
	if (g_hMmbCmdsPath == NULL)
		return FALSE;
	SendMessageW(g_hMmbCmdsPath, WM_SETFONT, (WPARAM)g_hMainFont, TRUE);

	g_hMmbCmdsPathBrowse = CreateWindowW(
		L"BUTTON",
		L"Browse...",
		WS_TABSTOP | WS_VISIBLE | WS_CHILD | BS_DEFPUSHBUTTON,
		rect.right - UI_BTN_WIDTH - UI_SIDE_OFFSET, y, UI_BTN_WIDTH, g_captionHighest,
		hWnd,
		(HMENU)IDC_BROWSE_CMDS_PATH,
		hInstance,
		NULL);
	if (g_hMmbCmdsPathBrowse == NULL)
		return FALSE;
	SendMessageW(g_hMmbCmdsPathBrowse, WM_SETFONT, (WPARAM)g_hMainFont, TRUE);

	// Output dir line
	y += g_formStep;

	g_hMmbWorkDirCap = CreateWindowW(
		L"STATIC",
		g_szMmbWorkDirCap,
		WS_VISIBLE | WS_CHILD | SS_CENTER,
		xL, y + (UI_HGT_PAD / 2) - 1, TextWidth(hWnd, g_szMmbWorkDirCap), g_captionHighest,
		hWnd,
		NULL,
		hInstance,
		NULL
	);
	if (g_hMmbWorkDirCap == NULL)
		return FALSE;
	SendMessageW(g_hMmbWorkDirCap, WM_SETFONT, (WPARAM)g_hMainFont, TRUE);

	g_hMmbWorkDir = CreateWindowW(
		L"EDIT",
		g_mmbWorkDir,
		WS_BORDER | WS_CHILD | WS_VISIBLE | ES_LEFT,
		xR, y, rect.right - xR - UI_CAPTION_SEP - UI_BTN_WIDTH - UI_SIDE_OFFSET, g_captionHighest,
		hWnd,
		(HMENU)ID_OUTPUT_DIR,
		hInstance,
		NULL
	);
	if (g_hMmbWorkDir == FALSE)
		return FALSE;
	SendMessageW(g_hMmbWorkDir, WM_SETFONT, (WPARAM)g_hMainFont, TRUE);

	g_hOutputDirBrowse = CreateWindowW(
		L"BUTTON",
		L"Browse...",
		WS_TABSTOP | WS_VISIBLE | WS_CHILD | BS_DEFPUSHBUTTON,
		rect.right - UI_BTN_WIDTH - UI_SIDE_OFFSET, y, UI_BTN_WIDTH, g_captionHighest,
		hWnd,
		(HMENU)IDC_BROWSE_OUTPUT_DIR,
		hInstance,
		NULL);
	if (g_hOutputDirBrowse == NULL)
		return FALSE;
	SendMessageW(g_hOutputDirBrowse, WM_SETFONT, (WPARAM)g_hMainFont, TRUE);

	// MMB output block
	y += g_formStep;

	g_hMmbOutputCap = CreateWindowW(
		L"STATIC",
		g_szMmbOutputCap,
		WS_VISIBLE | WS_CHILD,
		xL, y, g_captionWidest, g_captionHighest,
		hWnd,
		NULL,
		hInstance,
		NULL
	);
	if (g_hMmbOutputCap == NULL)
		return FALSE;
	SendMessageW(g_hMmbOutputCap, WM_SETFONT, (WPARAM)g_hMainFont, TRUE);

	y += g_formStep;
	LONG outputHeight = rect.bottom - 3 * (UI_BTN_HEIGHT - UI_CAPTION_SEP) - g_formStep - y - UI_CAPTION_SEP;

	g_hMmbOutput = CreateWindowW(
		L"EDIT",
		g_mmbOutput,
		WS_BORDER | WS_CHILD | WS_VISIBLE | WS_VSCROLL | WS_HSCROLL | ES_MULTILINE | ES_READONLY | ES_LEFT,
		xL, y, rect.right - 2 * UI_SIDE_OFFSET, outputHeight,
		hWnd,
		(HMENU)ID_MMB_OUTPUT,
		hInstance,
		NULL
	);
	SendMessageW(g_hMmbOutput, WM_SETFONT, (WPARAM)g_hMainFont, TRUE);

	y += outputHeight + UI_CAPTION_SEP;

	// Output controls line
	g_hDontScroll = CreateWindowW(
		L"BUTTON",
		L"Do not scroll",
		WS_VISIBLE | WS_CHILD | BS_CHECKBOX,
		xL, y, g_captionWidest, UI_BTN_HEIGHT,
		hWnd,
		(HMENU)IDC_DONT_SCROLL,
		hInstance,
		NULL
	);
	if (g_hDontScroll == NULL)
		return FALSE;
	CheckDlgButton(hWnd, IDC_DONT_SCROLL, g_scrollMmbOutput ? BST_UNCHECKED : BST_CHECKED);
	SendMessageW(g_hDontScroll, WM_SETFONT, (WPARAM)g_hMainFont, TRUE);

	g_hMmbOutputToFile = CreateWindowW(
		L"BUTTON",
		L"MMB output to file",
		WS_VISIBLE | WS_CHILD | BS_DEFPUSHBUTTON,
		xL + g_captionWidest + UI_CAPTION_SEP, y, g_captionWidest, UI_BTN_HEIGHT,
		hWnd,
		(HMENU)IDC_MMB_OUTPUT_TO_FILE,
		hInstance,
		NULL
	);
	if (g_hMmbOutputToFile == NULL)
		return FALSE;
	SendMessageW(g_hMmbOutputToFile, WM_SETFONT, (WPARAM)g_hMainFont, TRUE);

	y += UI_BTN_HEIGHT + UI_CAPTION_SEP;

	// MMB status line
	g_hMmbStatusCap = CreateWindowW(
		L"STATIC",
		g_szMmbStatusCap,
		WS_VISIBLE | WS_CHILD,
		xL, y, g_captionWidest, g_captionHighest,
		hWnd,
		NULL,
		hInstance,
		NULL
	);
	if (g_hMmbStatusCap == NULL)
		return FALSE;
	SendMessageW(g_hMmbStatusCap, WM_SETFONT, (WPARAM)g_hMainFont, TRUE);

	g_hMmbStatus = CreateWindowW(
		L"STATIC",
		MMB_ST_NOT_STARTED_CAP,
		WS_VISIBLE | WS_CHILD,
		xR, y, rect.right - xR, g_captionHighest,
		hWnd,
		NULL,
		hInstance,
		NULL
	);
	if (g_hMmbStatus == NULL)
		return FALSE;
	SendMessageW(g_hMmbStatus, WM_SETFONT, (WPARAM)g_hMainFont, TRUE);

	y = rect.bottom - UI_BTN_HEIGHT;
	LONG w = rect.right / 2 - UI_CAPTION_SEP - 2 * UI_SIDE_OFFSET;

	g_hRunBtn = CreateWindowW(
		L"BUTTON",
		L"Run",
		WS_TABSTOP | WS_VISIBLE | WS_CHILD | BS_DEFPUSHBUTTON,
		xL, y, xL + w, UI_BTN_HEIGHT,
		hWnd,
		(HMENU)IDC_RUN_MMB,
		hInstance,
		NULL
	);
	if (g_hRunBtn == NULL)
		return FALSE;
	SendMessageW(g_hRunBtn, WM_SETFONT, (WPARAM)g_hMainFont, TRUE);

	g_hStopBtn = CreateWindowW(
		L"BUTTON",
		L"Stop",
		WS_TABSTOP | WS_VISIBLE | WS_CHILD | BS_DEFPUSHBUTTON,
		rect.right - w - UI_CAPTION_SEP, y - UI_SIDE_OFFSET, w, UI_BTN_HEIGHT,
		hWnd,
		(HMENU)IDC_STOP_MMB,
		hInstance,
		NULL
	);
	if (g_hStopBtn == NULL)
		return FALSE;
	SendMessageW(g_hStopBtn, WM_SETFONT, (WPARAM)g_hMainFont, TRUE);

	return TRUE;
}

static
DWORD __stdcall MonitoringFunc(LPVOID param)
{
	DWORD exitCode;
	HWND hWnd = (HWND)param;

	assert(g_mmbProc.hProcess != NULL);

        BOOL run = TRUE;
	while (run) {
		char buf[512];

                run = WaitForSingleObject(g_mmbProc.hProcess, 0) == WAIT_TIMEOUT;

                DWORD read;
                do {
                        if (!ReadFile(g_hMmbReadPipe, buf, 512, &read, NULL))
                                break;

                        if (read > 0)
                                AppendMmbOutput(hWnd, buf, read);
                } while (read == 512);
	}

	GetExitCodeProcess(g_mmbProc.hProcess, &exitCode);

	CloseInvalidateHandle(g_hMmbReadPipe);

	PostMessageW(hWnd, WM_MMB_QUIT, exitCode, 0);

        return 0;
}

static
BOOL MonitorMmbProcess(HWND hWnd)
{
        assert(g_hMonitorThread == NULL);

	DWORD threadId;
	g_hMonitorThread = CreateThread(
		NULL,
		0,
		MonitoringFunc,
		(LPVOID)hWnd,
		0,
		&threadId
	);

	return g_hMonitorThread != NULL;
}

static
VOID ResizeControls(HWND hWnd)
{
	assert(g_formStep > 0);

	RECT rect;
	if (GetClientRect(hWnd, &rect) == FALSE)
		return;

	const LONG xL = UI_SIDE_OFFSET;
	const LONG xR = xL + g_captionWidest + UI_CAPTION_SEP;
	LONG y = UI_TOP_OFFSET;

	// MMB rubtime path line
	MoveWindow(
		g_hMmbRtDir,
		xR,  y, rect.right - xR - UI_BTN_WIDTH - UI_CAPTION_SEP - UI_SIDE_OFFSET, g_captionHighest,
		FALSE
	);

	MoveWindow(
		g_hMmbRtDirBrowse,
		rect.right - UI_BTN_WIDTH - UI_SIDE_OFFSET, y, UI_BTN_WIDTH, g_captionHighest,
		FALSE
	);

	// MMB commands line
	y += g_formStep;

	MoveWindow(
		g_hMmbCmdsPath,
		xR,  y, rect.right - xR - UI_BTN_WIDTH - UI_CAPTION_SEP - UI_SIDE_OFFSET, g_captionHighest,
		FALSE
	);

	MoveWindow(
		g_hMmbCmdsPathBrowse,
		rect.right - UI_BTN_WIDTH - UI_SIDE_OFFSET, y, UI_BTN_WIDTH, g_captionHighest,
		FALSE
	);

	// Output dir line
	y += g_formStep;

	MoveWindow(
		g_hMmbWorkDir,
		xR,  y, rect.right - xR - UI_BTN_WIDTH - UI_CAPTION_SEP - UI_SIDE_OFFSET, g_captionHighest,
		FALSE
	);

	MoveWindow(
		g_hOutputDirBrowse,
		rect.right - UI_BTN_WIDTH - UI_SIDE_OFFSET, y, UI_BTN_WIDTH, g_captionHighest,
		FALSE
	);

	// MMB output block
	y += 2 * g_formStep;
	LONG outputHeight = rect.bottom - 3 * (UI_BTN_HEIGHT - UI_CAPTION_SEP) - g_formStep - y - UI_CAPTION_SEP;

	MoveWindow(
		g_hMmbOutput,
		xL,  y, rect.right - 2 * UI_SIDE_OFFSET, outputHeight,
		FALSE
	);

	// Output controls line
	y += outputHeight + UI_CAPTION_SEP;


	MoveWindow(g_hDontScroll,
		xL, y, g_captionWidest, UI_BTN_HEIGHT,
		FALSE
	);

	MoveWindow(g_hMmbOutputToFile,
		xL + g_captionWidest + UI_CAPTION_SEP, y, g_captionWidest, UI_BTN_HEIGHT,
		FALSE
	);

	// MMB status line
	y += UI_BTN_HEIGHT + UI_CAPTION_SEP;

	MoveWindow(
		g_hMmbStatusCap,
		xL, y, g_captionWidest, g_captionHighest,
		FALSE
	);

	MoveWindow(
		g_hMmbStatus,
		xR, y, rect.right - xR, g_captionHighest,
		FALSE
	);

	// Run / stop buttons
	y = rect.bottom - UI_BTN_HEIGHT;
	LONG w = rect.right / 2 - UI_CAPTION_SEP - 2 * UI_SIDE_OFFSET;

	MoveWindow(
		g_hRunBtn,
		xL, y, xL + w, UI_BTN_HEIGHT,
		FALSE
	);

	MoveWindow(
		g_hStopBtn,
		rect.right - w - UI_CAPTION_SEP, y, w, UI_BTN_HEIGHT,
		FALSE
	);

	InvalidateRect(hWnd, &rect, TRUE);
}

static
VOID StartMmb(HWND hWnd)
{
	LPSTR mmbExec = NULL;
	LPSTR cmdLine = NULL;
	LPSTR mmbWorkDir = NULL;

	if (g_mmbProc.hProcess != NULL) {
		MessageBoxW(hWnd, L"MMB is already running", L"Cannot start MMB", MB_OK | MB_ICONINFORMATION);
		return;
	}

	if (wcslen(g_cmdsPath) == 0) {
		MessageBoxW(hWnd, L"Path to MMB commands file is empty", L"Cannot start MMB", MB_OK | MB_ICONEXCLAMATION);
		return;
	}

	if (!CheckMmbRtDirectory(g_mmbRtPath)) {
		MessageBoxW(
			hWnd,
			L"Path to the MMB runtime directory does not seem to contain usable MMB installation.",
			L"Input error",
			MB_OK | MB_ICONERROR
		);
		return;
	}

	if (wcslen(g_mmbWorkDir) < 1) {
		int answer = MessageBoxW(
			hWnd,
			L"MMB working directory is not set. This means that all files created by MMB will be "
			L"written into the MMB runtime directory. This is probably NOT what you want.\n\n"
			L"Do you want to continue?",
			L"Possibly wrong input",
			MB_YESNO | MB_ICONQUESTION
		);

		if (answer != IDYES)
			return;
	} else {
		DWORD error;
		if (!CopyMmbParametersFile(g_mmbRtPath, g_mmbWorkDir, &error)) {
			LPWSTR errorMsg = PrettyErrorString(
				L"Set MMB working directory does not contain \"parameter.csv\" file. Since this file must be present, "
				L"I tried to copy it there for you but the operation failed. Aborting...\n\n"
				L"Error",
				error
			);
			MessageBoxW(
				hWnd,
				errorMsg,
				L"I/O error",
				MB_OK | MB_ICONERROR
			);
			HFree(errorMsg, g_heap);
			return;
		}
	}

        mmbExec = GetMmbExecutable(g_mmbRtPath);
	if (mmbExec == NULL) {
		MessageBoxW(hWnd, L"Failed to get path to MMB executable", L"Cannot start MMB", MB_OK | MB_ICONERROR);
		goto out;
	}

	cmdLine = BuildMmbCommandLine(g_cmdsPath, g_mmbWorkDir);
	if (cmdLine == NULL) {
		MessageBoxW(hWnd, L"Failed to build MMB command line", L"Cannot start MMB", MB_OK | MB_ICONERROR);
		goto out;
	}

	mmbWorkDir = WideCharToSystemAnsi(g_mmbRtPath);
	if (mmbWorkDir == NULL) {
		MessageBoxW(hWnd, L"Failed to convert MMB working directory path", L"Cannot start MMB", MB_OK | MB_ICONERROR);
		goto out;
	}

	SECURITY_ATTRIBUTES secAttr = { sizeof(SECURITY_ATTRIBUTES) };
	secAttr.bInheritHandle = TRUE;
	secAttr.lpSecurityDescriptor = NULL;

        HANDLE hWrPipe;

        if (!CreatePipe(&g_hMmbReadPipe, &hWrPipe, &secAttr, 512)) {
		LPWSTR error = PrettyErrorString(L"Failed to create pipe", GetLastError());
		MessageBoxW(hWnd, error, L"Cannot start MMB", MB_OK | MB_ICONERROR);
		HFree(error, g_heap);

		goto out;
        }

	if (!SetHandleInformation(g_hMmbReadPipe, HANDLE_FLAG_INHERIT, 0)) {
		LPWSTR error = PrettyErrorString(L"Failed to set pipe handle information", GetLastError());
		MessageBoxW(hWnd, error, L"Cannot start MMB", MB_OK | MB_ICONERROR);
		HFree(error, g_heap);

		CloseInvalidateHandle(g_hMmbReadPipe);
                CloseHandle(hWrPipe);
		goto out;
	}

	ClearMmbOutput(hWnd);

	STARTUPINFOA si = { sizeof(STARTUPINFOA) };
        si.hStdError = hWrPipe;
        si.hStdOutput = hWrPipe;
	si.dwFlags |= STARTF_USESTDHANDLES;

	if (!CreateProcessA(
		mmbExec, cmdLine,
		NULL, NULL,
		TRUE,
		CREATE_NO_WINDOW,
		NULL,
		mmbWorkDir,
		&si,
		&g_mmbProc)) {
		LPWSTR error = PrettyErrorString(L"Failed to start MMB", GetLastError());
		MessageBoxW(hWnd, error, L"Cannot start MMB", MB_OK | MB_ICONERROR);
		HFree(error, g_heap);

                CloseInvalidateHandle(g_hMmbReadPipe);
                CloseHandle(hWrPipe);
		goto out;
	}

        CloseHandle(hWrPipe);

	SetMmbStatus(MMB_ST_RUNNING);

	MonitorMmbProcess(hWnd);

out:
	if (mmbWorkDir) HFree(mmbWorkDir, g_heap);
	if (cmdLine) HFree(cmdLine, g_heap);
	if (mmbExec) HFree(mmbExec, g_heap);
}

static
VOID StopMmb(HWND hWnd)
{
	if (g_mmbProc.hProcess == NULL)
		return;

	if (!TerminateProcess(g_mmbProc.hProcess, 1)) {
		LPWSTR error = PrettyErrorString(L"Failed to terminate MMB", GetLastError());
		MessageBoxW(hWnd, L"Failed to terminate MMB", L"Error", MB_OK | MB_ICONERROR);
		HFree(error, g_heap);
		return;
	}

	if (g_mmbProc.hProcess != NULL) {
		WaitForSingleObject(g_mmbProc.hProcess, INFINITE);
		CloseInvalidateHandle(g_mmbProc.hProcess);
		CloseInvalidateHandle(g_mmbProc.hThread);
	}
}

static
BOOL BrowseForFile(HWND hWnd, LPCWSTR initialPath, LPWSTR *outputPath, LPCWSTR title, BOOL readOnly)
{
        OPENFILENAMEW ofn = { sizeof(OPENFILENAMEW) };

        ofn.hwndOwner = hWnd;
	ofn.Flags |= OFN_PATHMUSTEXIST | OFN_NONETWORKBUTTON | (readOnly * OFN_FILEMUSTEXIST);
        ofn.nMaxFile = MAX_PATH;
        ofn.lpstrFile = HAlloc(WLEN(MAX_PATH), g_heap);
        ofn.lpstrTitle = title;
        ofn.lpstrFilter = L"Text\0*.txt;.dat\0All\0*.*\0\0";

        if (wcslen(initialPath) > 0)
            wcscpy_s(ofn.lpstrFile, ofn.nMaxFile, initialPath);
        else
            ofn.lpstrFile[0] = 0;

        BOOL ret = GetOpenFileNameW(&ofn);
        if (ret) {
		HFree(*outputPath, g_heap);

                DWORD len = wcslen(ofn.lpstrFile) + 1;
		*outputPath = HAlloc(WLEN(len), g_heap);
		wcscpy_s(*outputPath, len, ofn.lpstrFile);
	}

	HFree(ofn.lpstrFile, g_heap);

	return ret;
}

static
BOOL BrowseForFolder(HWND hWnd, LPWSTR *outputPath, LPCWSTR title)
{
	BROWSEINFOW bi = { hWnd };

	bi.hwndOwner = hWnd;
	bi.ulFlags = BIF_RETURNONLYFSDIRS | BIF_NEWDIALOGSTYLE;
	bi.lpszTitle = title;
	bi.iImage = -1;
	bi.pszDisplayName = HAlloc(WLEN(MAX_PATH), g_heap);
        bi.pszDisplayName[0] = 0;

	PIDLIST_ABSOLUTE pidl = SHBrowseForFolderW(&bi);

	if (pidl) {
		LPWSTR buf = HAlloc(WLEN(MAX_PATH), g_heap);
		BOOL ret = SHGetPathFromIDListW(pidl, buf);
		if (ret) {
			HFree(*outputPath, g_heap);
			*outputPath = buf;
		}

                HFree(bi.pszDisplayName, g_heap);
		CoTaskMemFree(pidl);
		return ret;
	}

	HFree(bi.pszDisplayName, g_heap);

	return FALSE;
}

static
BOOL MeasureCaptions(HWND hWnd)
{
        HDC hDc = NULL;

	hDc = GetDC(hWnd);
	if (hDc == NULL)
		return FALSE;

	SelectObject(hDc, g_hMainFont);

	for (size_t idx = 0; idx < ARRAY_SZ(g_formCaps); idx++) {
		LPCWSTR str = g_formCaps[idx];
		SIZE sz;
		if (GetTextExtentPointW(hDc, str, wcslen(str), &sz) == FALSE) {
			ReleaseDC(hWnd, hDc);
			return FALSE;
		}

		if (sz.cx > g_captionWidest)
			g_captionWidest = sz.cx;
		if (sz.cy > g_captionHighest)
			g_captionHighest = sz.cy;
	}

	ReleaseDC(hWnd, hDc);

	g_captionHighest += UI_HGT_PAD;

	g_formStep = g_captionHighest + 8;

	return TRUE;
}

static
VOID MmbOutputToFile(HWND hWnd)
{
	if (WaitForSingleObject(g_hOutputMutex, INFINITE) != WAIT_OBJECT_0)
		return;
	if (wcslen(g_mmbOutput) < 1) {
		ReleaseMutex(g_hOutputMutex);

		MessageBoxW(hWnd, L"MMB output is empty", L"No data", MB_OK | MB_ICONINFORMATION);
		return;
	}
	ReleaseMutex(g_hOutputMutex);

	LPWSTR path = NULL;

	if (!BrowseForFile(hWnd, L"", &path, L"Choose output file", FALSE))
		return;

	HANDLE hFile = CreateFileW(path, GENERIC_WRITE, 0, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
	if (hFile == INVALID_HANDLE_VALUE) {
		LPWSTR error = PrettyErrorString(L"Failed to open output file", GetLastError());
		MessageBoxW(hWnd, L"Failed to open output file", L"I/O error", MB_OK | MB_ICONERROR);
		HFree(error, g_heap);

		HFree(path, g_heap);
		return;
	}
	HFree(path, g_heap);

	if (WaitForSingleObject(g_hOutputMutex, INFINITE) != WAIT_OBJECT_0)
		return;

	size_t len = WLEN(wcslen(g_mmbOutput));
	DWORD written = 0;

	while (written < len) {
		if (!WriteFile(hFile, g_mmbOutput + written, len - written, &written, NULL)) {
			LPWSTR error = PrettyErrorString(L"Failed to write output", GetLastError());
			MessageBoxW(hWnd, L"Failed to write output", L"I/O error", MB_OK | MB_ICONERROR);
			HFree(error, g_heap);

			break;
		}
	}

	CloseHandle(hFile);
	ReleaseMutex(g_hOutputMutex);
}

static
VOID ToggleMmbOutputScroll(HWND hWnd)
{
	BOOL checked = IsDlgButtonChecked(hWnd, IDC_DONT_SCROLL);

	CheckDlgButton(hWnd, IDC_DONT_SCROLL, !checked ? BST_CHECKED : BST_UNCHECKED);
	g_scrollMmbOutput = checked;
}

static
VOID UpdateFormInput(HWND hWnd, LPWSTR *outputStr)
{
	DWORD len = GetWindowTextLengthW(hWnd);

        HFree(*outputStr, g_heap);
	*outputStr = HAlloc(WLEN(len + 1), g_heap);
	GetWindowTextW(hWnd, *outputStr, len + 1);
}

static
VOID UpdateMmbOutput()
{
	if (WaitForSingleObject(g_hOutputMutex, INFINITE) != WAIT_OBJECT_0)
		return;

	SetWindowTextW(g_hMmbOutput, g_mmbOutput);

	ReleaseMutex(g_hOutputMutex);

	if (g_scrollMmbOutput)
		PostMessageW(g_hMmbOutput, WM_VSCROLL, LOWORD(SB_BOTTOM), 0);
}

static
LRESULT CALLBACK WndProc(HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam)
{
	switch (msg) {
	case WM_CREATE:
		SendMessageW(hWnd, WM_SETFONT, (WPARAM)g_hMainFont, TRUE);

		if (MeasureCaptions(hWnd) == FALSE) {
			MessageBoxW(hWnd, L"Cannot calculate controls layout", L"Initialization error", MB_OK | MB_ICONERROR);
			PostQuitMessage(0);
			return 0;
		}

		if (InitControls(hWnd, ((LPCREATESTRUCTW)lParam)->hInstance) == FALSE) {
			MessageBoxW(hWnd, L"Cannot create controls", L"Initialization error", MB_OK | MB_ICONERROR);
			PostQuitMessage(0);
			return 0;
		}

		break;
	case WM_COMMAND:
	{
		const DWORD hi = HIWORD(wParam);
		const DWORD lo = LOWORD(wParam);

		switch (hi) {
		case 1:
			switch (lo) {
			case IDA_RUN_MMB:
				StartMmb(hWnd);
				break;
			case IDA_STOP_MMB:
				StopMmb(hWnd);
				break;
			}
			break;
		case EN_CHANGE:
			switch (lo) {
			case ID_MMB_RT:
				UpdateFormInput((HWND)lParam, &g_mmbRtPath);
				break;
			case ID_MMB_CMDS_PATH:
				UpdateFormInput((HWND)lParam, &g_cmdsPath);
				break;
			case ID_OUTPUT_DIR:
				UpdateFormInput((HWND)lParam, &g_mmbWorkDir);
				break;
			}
			break;
		case BN_CLICKED:
			switch (lo) {
			case IDC_BROWSE_MMB_RT:
				if (BrowseForFolder(hWnd, &g_mmbRtPath, L"Choose MMB runtime folder"))
					SetWindowTextW(g_hMmbRtDir, g_mmbRtPath);
				break;
			case IDC_BROWSE_CMDS_PATH:
				if (BrowseForFile(hWnd, g_cmdsPath, &g_cmdsPath, L"Choose MMB commands file", TRUE))
					SetWindowTextW(g_hMmbCmdsPath, g_cmdsPath);
				break;
			case IDC_BROWSE_OUTPUT_DIR:
				if (BrowseForFolder(hWnd, &g_mmbWorkDir, L"Choose MMB working directory"))
					SetWindowTextW(g_hMmbWorkDir, g_mmbWorkDir);
				break;
			case IDC_RUN_MMB:
				StartMmb(hWnd);
				break;
			case IDC_STOP_MMB:
				StopMmb(hWnd);
				break;
			case IDC_DONT_SCROLL:
				ToggleMmbOutputScroll(hWnd);
				break;
			case IDC_MMB_OUTPUT_TO_FILE:
				MmbOutputToFile(hWnd);
				break;
			}
			break;
		}
		break;
	}
	case WM_SIZE:
		ResizeControls(hWnd);
		break;
	case WM_DESTROY:
		StopMmb(hWnd);
		PostQuitMessage(0);
		break;
	case WM_MMB_QUIT:
		FinalizeMmb(hWnd, wParam);
		break;
	case WM_MMB_OUTPUT:
		UpdateMmbOutput();
		break;
	default:
		return DefWindowProcW(hWnd, msg, wParam, lParam);
	}

	return 0;
}

static
BOOL RegisterWinMmbWindowClass(HINSTANCE hInstance)
{
	WNDCLASSEXW wcex = { sizeof(WNDCLASSEXW) };
	wcex.style = 0;
	wcex.lpfnWndProc = WndProc;
	wcex.cbClsExtra = 0;
	wcex.cbWndExtra = 0;
	wcex.hInstance = hInstance;
	wcex.hIcon = LoadIconW(hInstance, MAKEINTRESOURCEW(IDI_ICON1));
	wcex.hCursor = LoadCursorW(NULL, MAKEINTRESOURCEW(32512));
	wcex.hbrBackground = (HBRUSH)(COLOR_WINDOW);
	wcex.lpszMenuName = NULL;
	wcex.lpszClassName = g_szWindowClass;
	wcex.hIconSm = LoadIconW(wcex.hInstance, MAKEINTRESOURCEW(IDI_ICON1));

	return RegisterClassExW(&wcex);
}

#ifdef _MSC_VER
int APIENTRY wWinMain(
	_In_ HINSTANCE hInstance,
	_In_opt_ HINSTANCE hPrevInstance,
	_In_ LPWSTR lpCmdLine,
	_In_ int nCmdShow)
#else
int APIENTRY WinMain(
	HINSTANCE hInstance,
	HINSTANCE hPrevInstance,
	LPSTR     lpCmdLine,
	int       nCmdShow
)
#endif // __MSC_VER
{
	UNREFERENCED_PARAMETER(hPrevInstance);
	UNREFERENCED_PARAMETER(lpCmdLine);

	INITCOMMONCONTROLSEX icc;
	icc.dwSize = sizeof(icc);
	icc.dwICC = ICC_WIN95_CLASSES;
	InitCommonControlsEx(&icc);

	g_heap = GetProcessHeap();
	if (g_heap == NULL)
		return EXIT_FAILURE;

	HACCEL hAccelTable = LoadAcceleratorsW(hInstance, L"Accels");
	if (hAccelTable == NULL)
		return EXIT_FAILURE;

	if (!InitMainFont())
		return EXIT_FAILURE;

	ZeroMemory(&g_mmbProc, sizeof(PROCESS_INFORMATION));

        g_hOutputMutex = CreateMutexW(NULL, FALSE, L"MmbOutputMutex");
        if (g_hOutputMutex == NULL) {
            MessageBoxW(NULL, L"Failed to create mutex", L"Initialization error", MB_OK | MB_ICONERROR);
            return EXIT_FAILURE;
        }

	if (!InitializeDefaults()) {
		MessageBoxW(NULL, L"Failed to initialize application defaults", L"Initialization error", MB_OK | MB_ICONERROR);
		return EXIT_FAILURE;
	}

	if (!RegisterWinMmbWindowClass(hInstance)) {
		MessageBoxW(NULL, L"Failed to register window class", L"Initialization error", MB_OK | MB_ICONERROR);

		return EXIT_FAILURE;
	}

	HWND hWnd = CreateWindowW(
		g_szWindowClass,
		g_szWindowTitle,
		WS_OVERLAPPEDWINDOW,
		CW_USEDEFAULT,
		CW_USEDEFAULT,
		CW_USEDEFAULT,
		CW_USEDEFAULT,
		NULL,
		NULL,
		hInstance,
		NULL
	);
	if (hWnd == NULL) {
		MessageBoxW(NULL, L"Failed to create window", L"Initialization error", MB_OK | MB_ICONERROR);
		return EXIT_FAILURE;
	}

	ShowWindow(hWnd, nCmdShow);
	UpdateWindow(hWnd);

	MSG msg;
	while (GetMessageW(&msg, NULL, 0, 0)) {
		if (!TranslateAcceleratorW(msg.hwnd, hAccelTable, &msg)) {
			TranslateMessage(&msg);
			DispatchMessageW(&msg);
		}
	}

	FreeAllocated();

	return msg.wParam;
}
