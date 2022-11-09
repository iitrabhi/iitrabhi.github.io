---
layout: post
title: "How to get 10 GB of sync space for your zotero library."
description: "Pcloud ended support for webdav sync in Zotero. It's time to move to box.com"
categories: [productivity, presentation]
typora-root-url: ../../../../website
---

## The problem

For the past three years I was using Pcloud to sync my zotero library across devices. It provided 10GB of free space that was more than enough for all of my documents. But since feburary 2022, pcloud ended support for webdav sync and thus now we can not sync our library between office and home system.

## The solution - Box.com

Box.com is a service similar to pcloud that offers 10GB of free space. Follow these steps to setup box.com for zotero on MAC, Windows and Ipad

- Make an account on [Box.com](https://www.box.com)

- Goto Zotero →  Preferences → Sync

- Under file syncing choose WebDAV

- URL - `dav.box.com/dav`

- Enter your credentials for Box.com in username and password

  ![image-20220428180050793](/assets/images/image-20220428180050793.png)

- For windows and ipad follow the same process. Now, you can enjoy free sync of your complete library between devices.

## PS

If after changing your WEBDAV server your files are not syncing then goto Zotero →  Preferences → Sync → Reset. Choose the radio button Replace Online Library and hit Reset.

`Caution: The above process will overwrite your online library with your local library. Use this option wisely.`

## References

-  [Pcloud end webdav support ](https://forums.zotero.org/discussion/94430/pcloud-ending-webdav-support-for-free-plan)
- [Zotero forum box setup](https://forums.zotero.org/discussion/33934/file-sync-over-webdav-to-box-com)

